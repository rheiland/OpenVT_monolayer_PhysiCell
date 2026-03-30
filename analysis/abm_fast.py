"""
Growing monolayer model in Python  — NumPy-accelerated version
==============================================================

Optimisations over abm_slow.py
-------------------------------
1. Pre-allocated numpy arrays of shape (MAX_AGENTS,) for every per-agent
   field (x, y, vel_x, vel_y, prev_vel_x, prev_vel_y, area, growth_rate,
   division_area, norm_rand, ID).  Only the first `n` rows are "live".
2. Growth step: a single vectorised `area[0:n] += growth_rate[0:n]`.
3. Division step: a boolean mask finds all dividers; daughters are written
   into a contiguous block at the end of the live region — no Python-level
   list allocation.
4. Force / velocity step: fully vectorised pairwise computation using
   numpy broadcasting (N×N distance matrix).  For very large N an optional
   spatial-hash neighbour list can be dropped in, but broadcasting is already
   ~100× faster than nested Python loops for N ≤ a few thousand.
5. Frame storage: each frame is saved as a dict of numpy slices (views /
   lightweight copies) rather than deepcopy-ing a list of dataclass objects.
"""

import sys
import math
import random
import csv

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

print("# args=", len(sys.argv))
if len(sys.argv) < 5:
    print("Usage: <repulsion(10)> <max_cells> <win_size> <split_type(0-2)>\n")
    sys.exit(1)

idx = 1
repulsion      = float(sys.argv[idx]); idx += 1
MAX_AGENTS     = int(sys.argv[idx]);   idx += 1
win_size       = float(sys.argv[idx]); idx += 1
daughter_split = int(sys.argv[idx])

# ---------------------------------------------------------------------------
# SimState  — thin wrapper around pre-allocated arrays
# ---------------------------------------------------------------------------

class SimState:
    """
    All per-agent data stored in flat numpy arrays.
    Only indices [0 : self.n] are considered live agents.
    """

    def __init__(self, max_agents: int, n_seed: int,
                 growth_rate: float, division_area: float,
                 norm_rand: float, initial_area: float):

        cap = max_agents * 2 + 16   # a little headroom for one step's worth of divisions

        self.cap = cap
        self.n   = 0   # live count

        self.x          = np.zeros(cap, dtype=np.float64)
        self.y          = np.zeros(cap, dtype=np.float64)
        self.vel_x      = np.zeros(cap, dtype=np.float64)
        self.vel_y      = np.zeros(cap, dtype=np.float64)
        self.prev_vel_x = np.zeros(cap, dtype=np.float64)
        self.prev_vel_y = np.zeros(cap, dtype=np.float64)
        self.area        = np.zeros(cap, dtype=np.float64)
        self.growth_rate = np.zeros(cap, dtype=np.float64)
        self.div_area    = np.zeros(cap, dtype=np.float64)
        self.norm_rand   = np.zeros(cap, dtype=np.float64)
        self.agent_id    = np.zeros(cap, dtype=np.int64)

        self._next_id = 0
        self.time     = 0

        # seed agents
        for _ in range(n_seed):
            self._add_agent(x=0.0, y=0.0, area=initial_area,
                            growth_rate=growth_rate,
                            division_area=division_area,
                            norm_rand=norm_rand)

    # ------------------------------------------------------------------
    def _add_agent(self, x, y, area, growth_rate, division_area, norm_rand,
                   reuse_idx: int = -1):
        """Write a new agent into slot reuse_idx (for daughter-1) or append."""
        if reuse_idx >= 0:
            i = reuse_idx
        else:
            i = self.n
            self.n += 1

        self.x[i]           = x
        self.y[i]           = y
        self.vel_x[i]       = 0.0
        self.vel_y[i]       = 0.0
        self.prev_vel_x[i]  = 0.0
        self.prev_vel_y[i]  = 0.0
        self.area[i]        = area
        self.growth_rate[i] = growth_rate
        self.div_area[i]    = division_area
        self.norm_rand[i]   = norm_rand
        self.agent_id[i]    = self._next_id
        self._next_id      += 1

    # ------------------------------------------------------------------
    def _snapshot(self) -> dict:
        """Lightweight copy of the live region — used for frame storage."""
        n = self.n
        return dict(
            x        = self.x[:n].copy(),
            y        = self.y[:n].copy(),
            area     = self.area[:n].copy(),
            norm_rand= self.norm_rand[:n].copy(),
        )

    # ------------------------------------------------------------------
    def _draw_rand_positive(self) -> float:
        v = random.normalvariate(2.0, 0.4)
        while v <= 0:
            v = random.normalvariate(2.0, 0.4)
        return v

    # ------------------------------------------------------------------
    def step(self) -> int:
        n = self.n

        # 1) Vectorised growth
        self.area[:n] += self.growth_rate[:n]

        # 2) Find dividers
        dividers = np.where(self.area[:n] >= self.div_area[:n])[0]

        # Process divisions (in reverse index order so we can safely append)
        divisions = 0
        for i in dividers:
            da = self.area[i] / 2.0
            dr = math.sqrt(da / math.pi)
            theta = random.random() * 6.283185307179586
            xvec  = math.cos(theta) * dr * 2
            yvec  = math.sin(theta) * dr * 2

            nr1 = self._draw_rand_positive()
            nr2 = self._draw_rand_positive()
            div1 = 78.54 * nr1
            div2 = 78.54 * nr2

            gr = self.growth_rate[i]

            if daughter_split == 0:
                dx1, dy1 = self.x[i],          self.y[i]
                dx2, dy2 = self.x[i] + xvec,   self.y[i] + yvec
            elif daughter_split == 1:
                dx1, dy1 = self.x[i] + xvec/2, self.y[i] + yvec/2
                dx2, dy2 = self.x[i] - xvec/2, self.y[i] - yvec/2
            else:   # == 2
                dx1, dy1 = self.x[i] + xvec/2.5, self.y[i] + yvec/2.5
                dx2, dy2 = self.x[i] - xvec/2.5, self.y[i] - yvec/2.5

            # Daughter-1 reuses slot i; daughter-2 is appended
            old_id = self.agent_id[i]
            self._add_agent(dx1, dy1, da, gr, div1, nr1, reuse_idx=i)
            self.agent_id[i] = old_id          # keep parent ID on d1
            self._add_agent(dx2, dy2, da, gr, div2, nr2)

            divisions += 1

        # 3) Vectorised relaxation (force / velocity / position)
        dt_mechanics         = 0.1
        max_relax_steps      = 10
        # argh: changing these doesn't seem to affect overcrowding, so...
        # dt_mechanics         = 0.05
        # max_relax_steps      = 20
        # dt_mechanics         = 0.02
        # max_relax_steps      = 50
        cc_repulsion         = repulsion
        d1_coef              = dt_mechanics * 1.5
        d2_coef              = dt_mechanics * -0.5
        n = self.n   # may have grown via divisions

        for _ in range(max_relax_steps):
            # --- compute pairwise displacement & distance ---
            xi = self.x[:n]
            yi = self.y[:n]

            dx  = xi[:, None] - xi[None, :]   # (n,n)
            dy  = yi[:, None] - yi[None, :]
            dist2 = dx*dx + dy*dy
            np.fill_diagonal(dist2, np.inf)    # exclude self-pairs
            dist  = np.sqrt(dist2)

            # combined radii  R_ij = r_i + r_j
            ri  = np.sqrt(self.area[:n] / math.pi)
            R   = ri[:, None] + ri[None, :]    # (n,n)

            # repulsion only when dist < 1.25*R
            # in_range  = dist < 1.25 * R
            in_range  = dist < 1.05 * R
            overlapping = dist < R

            # temp_r = (1 - d/R)^2 * repulsion   (only where overlapping)
            with np.errstate(invalid='ignore', divide='ignore'):
                ratio = np.where(overlapping & in_range,
                                 1.0 - dist / R, 0.0)
            temp_r = ratio * ratio * cc_repulsion   # (n,n)

            # force direction: normalise displacement
            with np.errstate(invalid='ignore', divide='ignore'):
                inv_dist = np.where(in_range & overlapping, 1.0 / dist, 0.0)

            fx = (dx * inv_dist * temp_r).sum(axis=1)   # (n,)
            fy = (dy * inv_dist * temp_r).sum(axis=1)

            # Adams-Bashforth position update
            self.x[:n]  += fx * d1_coef + self.prev_vel_x[:n] * d2_coef
            self.y[:n]  += fy * d1_coef + self.prev_vel_y[:n] * d2_coef

            self.prev_vel_x[:n] = fx
            self.prev_vel_y[:n] = fy
            self.vel_x[:n] = 0.0
            self.vel_y[:n] = 0.0

        self.time += 1
        return divisions


# ---------------------------------------------------------------------------
# Simulation — orchestrates SimState, stores frames
# ---------------------------------------------------------------------------

class Simulation:

    def __init__(self, n_seed=1, growth_rate=0.0, division_area=0.0,
                 norm_rand=0.0, initial_area=0.0):
        self._state = SimState(
            max_agents    = MAX_AGENTS,
            n_seed        = n_seed,
            growth_rate   = growth_rate,
            division_area = division_area,
            norm_rand     = norm_rand,
            initial_area  = initial_area,
        )

    @property
    def n_agents(self):
        return self._state.n

    def precompute(self, steps: int, max_agents: int = 500) -> None:
        import copy
        print(f"Pre-computing {steps} steps...", end=" ", flush=True)
        s = self._state
        self._frames: list[dict] = [s._snapshot()]

        for _ in range(steps):
            s.step()
            self._frames.append(s._snapshot())
            if s.n >= max_agents:
                print(f"(stopped early at t={s.time} — {s.n} agents)", end=" ")
                break

        print(f"done. {len(self._frames)} frames stored.")

        # Save final CSV
        # file_out = 'py_monolayer_small.csv'
        # print("--> ", file_out)
        # snap = self._frames[-1]
        # with open(file_out, "w", newline="") as fh:
        #     writer = csv.writer(fh)
        #     writer.writerow(['x_pos', 'y_pos', 'radius_i', 'norm_rand_i'])
        #     for xi, yi, ai, nri in zip(snap['x'], snap['y'],
        #                                snap['area'], snap['norm_rand']):
        #         writer.writerow([xi, yi, math.sqrt(ai / math.pi), nri])

    def interactive_viewer(self, title: str = "Agent-Based Model",
                           interval: int = 100):
        from matplotlib.widgets import Button, Slider
        from matplotlib.animation import FuncAnimation

        if not hasattr(self, "_frames"):
            raise RuntimeError("Call precompute() before interactive_viewer().")

        frames   = self._frames
        n_frames = len(frames)

        fig = plt.figure(figsize=(6, 6), facecolor="#0d0d0d")
        fig.suptitle(title, fontsize=13, fontweight="bold", color="white")

        ax = fig.add_axes([0.1, 0.2, 0.8, 0.7])
        ax.set_facecolor("#111111")
        ax.tick_params(colors="white")
        for sp in ax.spines.values():
            sp.set_edgecolor("#444444")
        ax.set_aspect("equal")
        ax.set_xlabel("x", color="white")
        ax.set_ylabel("y", color="white")
        ax.set_xlim(-win_size, win_size)
        ax.set_ylim(-win_size, win_size)
        title_text = ax.set_title("", color="white", fontsize=10)

        btn_color = "#1e2d3d"
        btn_hover = "#2e4d6d"

        ax_slider = fig.add_axes([0.10, 0.10, 0.80, 0.03], facecolor="#1a1a2e")
        ax_back   = fig.add_axes([0.28, 0.02, 0.08, 0.055], facecolor=btn_color)
        ax_play   = fig.add_axes([0.38, 0.02, 0.10, 0.055], facecolor=btn_color)
        ax_fwd    = fig.add_axes([0.50, 0.02, 0.08, 0.055], facecolor=btn_color)

        slider   = Slider(ax_slider, "t", 0, n_frames - 1,
                          valinit=0, valstep=1, color="#00d4ff")
        slider.label.set_color("white")
        slider.valtext.set_color("white")

        btn_back = Button(ax_back, "◀◀", color=btn_color, hovercolor=btn_hover)
        btn_play = Button(ax_play, "▶ Play", color=btn_color, hovercolor=btn_hover)
        btn_fwd  = Button(ax_fwd,  "▶▶", color=btn_color, hovercolor=btn_hover)
        for btn in (btn_back, btn_play, btn_fwd):
            btn.label.set_color("white")
            btn.label.set_fontsize(10)

        state = {"frame": 0, "playing": False}

        def _draw_frame(fidx: int) -> None:
            fidx = max(0, min(fidx, n_frames - 1))
            state["frame"] = fidx
            snap = frames[fidx]

            for p in list(ax.patches):
                p.remove()
            radii = np.sqrt(snap['area'] / math.pi)
            for xi, yi, ri in zip(snap['x'], snap['y'], radii):
                ax.add_patch(patches.Circle(
                    (xi, yi), radius=ri,
                    facecolor='gray', alpha=0.72,
                    linewidth=1.0, edgecolor="white",
                ))

            title_text.set_text(
                f"t = {fidx}   |   agents = {len(snap['x'])}"
            )
            slider.eventson = False
            slider.set_val(fidx)
            slider.eventson = True
            fig.canvas.draw_idle()

        def on_back(_e):
            state["playing"] = False; btn_play.label.set_text("▶ Play")
            _draw_frame(state["frame"] - 1)

        def on_fwd(_e):
            state["playing"] = False; btn_play.label.set_text("▶ Play")
            _draw_frame(state["frame"] + 1)

        def on_play(_e):
            state["playing"] = not state["playing"]
            btn_play.label.set_text("■ Pause" if state["playing"] else "▶ Play")
            fig.canvas.draw_idle()

        def on_slider(val):
            state["playing"] = False; btn_play.label.set_text("▶ Play")
            _draw_frame(int(val))

        btn_back.on_clicked(on_back)
        btn_fwd.on_clicked(on_fwd)
        btn_play.on_clicked(on_play)
        slider.on_changed(on_slider)

        def _tick(_f):
            if state["playing"]:
                nxt = state["frame"] + 1
                if nxt >= n_frames:
                    state["playing"] = False
                    btn_play.label.set_text("▶ Play")
                else:
                    _draw_frame(nxt)

        anim = FuncAnimation(fig, _tick, interval=interval,
                             cache_frame_data=False)
        _draw_frame(0)
        plt.show()
        return anim


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    random.seed(42)

    sim = Simulation(
        n_seed        = 1,
        growth_rate   = 0.1778,
        division_area = 157.08,
        norm_rand     = 2.0,
        initial_area  = 78.54,
    )

    sim.precompute(steps=50000, max_agents=MAX_AGENTS)
    sim.interactive_viewer(title="Growing Monolayer (vectorized)", interval=1)
