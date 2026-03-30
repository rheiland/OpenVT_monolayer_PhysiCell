"""
Chaste-Style Off-Lattice Cell Growth & Division Simulation
===========================================================
Implements the core Chaste abstractions:
  - CellPopulation  (NodeBasedCellPopulation equivalent)
  - Cell            (with a StochasticCellCycleModel)
  - Meineke linear spring force
  - OffLatticeSimulation loop

Starting condition: single cell at the origin.
Output: animated visualisation + summary plot.

Dependencies: numpy, matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import matplotlib.gridspec as gridspec
from dataclasses import dataclass, field
from typing import List, Optional
import random

# ──────────────────────────────────────────────────────────────────────────────
# PARAMETERS  (mirror Chaste defaults where possible)
# ──────────────────────────────────────────────────────────────────────────────
CELL_RADIUS          = 0.5          # equilibrium radius (Chaste: 0.5)
SPRING_STIFFNESS     = 15.0         # μ in Meineke model
DAMPING              = 1.0          # η (over-damped mechanics)
CUTOFF_MULTIPLE      = 1.5          # neighbour search radius = CUTOFF_MULTIPLE * 2R
DT                   = 0.005        # time-step (hours)
END_TIME             = 24.0         # total simulation time (hours)
MEAN_CELL_CYCLE_TIME = 4.0          # hours (stochastic Gaussian)
SD_CELL_CYCLE_TIME   = 0.5          # hours
BIRTH_SEPARATION     = 0.1          # daughter cells placed ±this from parent
SEED                 = 42
RECORD_EVERY         = 0.1          # hours between snapshots

np.random.seed(SEED)
random.seed(SEED)

# ──────────────────────────────────────────────────────────────────────────────
# CELL CYCLE MODEL
# ──────────────────────────────────────────────────────────────────────────────
class StochasticCellCycleModel:
    """
    Mirrors Chaste's AbstractSimpleCellCycleModel.
    Each cell draws a cycle duration from N(mean, sd²) at birth.
    ReadyToDivide() returns True once age >= cycle_duration.
    """
    def __init__(self, birth_time: float = 0.0):
        self.birth_time    = birth_time
        self.cycle_duration = max(
            2.0,
            np.random.normal(MEAN_CELL_CYCLE_TIME, SD_CELL_CYCLE_TIME)
        )

    @property
    def age(self) -> float:
        return _current_time - self.birth_time          # global clock

    def ready_to_divide(self) -> bool:
        return self.age >= self.cycle_duration

    def create_daughter_model(self, birth_time: float) -> "StochasticCellCycleModel":
        m = StochasticCellCycleModel(birth_time)
        return m


# ──────────────────────────────────────────────────────────────────────────────
# CELL
# ──────────────────────────────────────────────────────────────────────────────
_cell_id_counter = 0

@dataclass
class Cell:
    """Mirrors Chaste's CellPtr."""
    position     : np.ndarray
    cell_cycle   : StochasticCellCycleModel
    cell_id      : int   = field(default=None)
    generation   : int   = 0
    radius       : float = CELL_RADIUS

    def __post_init__(self):
        global _cell_id_counter
        if self.cell_id is None:
            self.cell_id = _cell_id_counter
            _cell_id_counter += 1

    def divide(self, current_time: float):
        """
        Returns two daughter Cells placed symmetrically around parent.
        Mirrors Chaste's CellDivisionRule.
        """
        angle    = np.random.uniform(0, 2 * np.pi)
        disp     = np.array([np.cos(angle), np.sin(angle)]) * BIRTH_SEPARATION
        d1 = Cell(
            position   = self.position + disp,
            cell_cycle = self.cell_cycle.create_daughter_model(current_time),
            generation = self.generation + 1,
        )
        d2 = Cell(
            position   = self.position - disp,
            cell_cycle = self.cell_cycle.create_daughter_model(current_time),
            generation = self.generation + 1,
        )
        return d1, d2


# ──────────────────────────────────────────────────────────────────────────────
# FORCE MODEL  –  Meineke Linear Spring
# ──────────────────────────────────────────────────────────────────────────────
def compute_forces(cells: List[Cell]) -> List[np.ndarray]:
    """
    Implements Chaste's GeneralisedLinearSpringForce.
    F_ij = μ (|r_ij| - s_ij) r̂_ij   if |r_ij| < cutoff
    where s_ij = R_i + R_j  (natural length = sum of radii).
    """
    n      = len(cells)
    forces = [np.zeros(2) for _ in range(n)]
    cutoff = CUTOFF_MULTIPLE * 2 * CELL_RADIUS

    for i in range(n):
        for j in range(i + 1, n):
            delta = cells[j].position - cells[i].position
            dist  = np.linalg.norm(delta)
            if dist == 0 or dist > cutoff:
                continue
            natural_len = cells[i].radius + cells[j].radius
            overlap     = dist - natural_len          # negative → compression
            direction   = delta / dist
            f           = SPRING_STIFFNESS * overlap * direction
            forces[i]  += f
            forces[j]  -= f
    return forces


# ──────────────────────────────────────────────────────────────────────────────
# OFF-LATTICE SIMULATION
# ──────────────────────────────────────────────────────────────────────────────
_current_time: float = 0.0

def run_simulation():
    global _current_time, _cell_id_counter
    _current_time    = 0.0
    _cell_id_counter = 0

    # Initialise with a single cell at the origin
    initial_cell = Cell(
        position   = np.array([0.0, 0.0]),
        cell_cycle = StochasticCellCycleModel(birth_time=0.0),
    )
    cells: List[Cell] = [initial_cell]

    snapshots        = []    # list of (time, positions, generations)
    cell_count_ts    = []    # (time, count)
    next_record_time = 0.0

    n_steps = int(END_TIME / DT)

    for step in range(n_steps):
        _current_time = step * DT

        # ── 1. Identify cells ready to divide ──
        to_divide = [c for c in cells if c.cell_cycle.ready_to_divide()]
        new_cells  = []
        for parent in to_divide:
            d1, d2 = parent.divide(_current_time)
            new_cells.extend([d1, d2])
        # Remove parents, add daughters
        if to_divide:
            divide_ids = {c.cell_id for c in to_divide}
            cells = [c for c in cells if c.cell_id not in divide_ids] + new_cells

        # ── 2. Compute forces (Meineke spring) ──
        forces = compute_forces(cells)

        # ── 3. Update positions (over-damped Euler: η dx/dt = F) ──
        for cell, force in zip(cells, forces):
            cell.position += (DT / DAMPING) * force

        # ── 4. Record snapshot ──
        if _current_time >= next_record_time:
            positions   = np.array([c.position   for c in cells])
            generations = np.array([c.generation for c in cells])
            ages        = np.array([c.cell_cycle.age for c in cells])
            snapshots.append((_current_time, positions.copy(),
                              generations.copy(), ages.copy()))
            cell_count_ts.append((_current_time, len(cells)))
            next_record_time += RECORD_EVERY

    return snapshots, cell_count_ts


# ──────────────────────────────────────────────────────────────────────────────
# VISUALISATION
# ──────────────────────────────────────────────────────────────────────────────
def visualise(snapshots, cell_count_ts):
    fig = plt.figure(figsize=(14, 8), facecolor="#0d1117")
    gs  = gridspec.GridSpec(2, 2, figure=fig,
                            left=0.06, right=0.97,
                            top=0.93,  bottom=0.08,
                            hspace=0.42, wspace=0.35)

    ax_main  = fig.add_subplot(gs[:, 0])   # large colony view
    ax_count = fig.add_subplot(gs[0, 1])   # cell count over time
    ax_gen   = fig.add_subplot(gs[1, 1])   # generation histogram

    for ax in [ax_main, ax_count, ax_gen]:
        ax.set_facecolor("#161b22")
        for spine in ax.spines.values():
            spine.set_edgecolor("#30363d")
        ax.tick_params(colors="#8b949e")
        ax.xaxis.label.set_color("#8b949e")
        ax.yaxis.label.set_color("#8b949e")
        ax.title.set_color("#e6edf3")

    # colour map: generation → colour
    max_gen = max(g.max() for _, _, g, _ in snapshots if len(g))
    cmap    = plt.cm.plasma
    norm    = Normalize(vmin=0, vmax=max(max_gen, 1))

    # ── Static plots ──
    times  = [t for t, _ in cell_count_ts]
    counts = [c for _, c in cell_count_ts]
    ax_count.plot(times, counts, color="#58a6ff", linewidth=1.8)
    ax_count.set_xlabel("Time (hours)")
    ax_count.set_ylabel("Cell count")
    ax_count.set_title("Cell Population Growth")
    ax_count.grid(True, color="#21262d", linewidth=0.5)

    final_gens = snapshots[-1][2]
    unique_gens, gen_counts = np.unique(final_gens, return_counts=True)
    bars = ax_gen.bar(unique_gens, gen_counts,
                      color=[cmap(norm(g)) for g in unique_gens],
                      edgecolor="#0d1117", linewidth=0.5)
    ax_gen.set_xlabel("Generation")
    ax_gen.set_ylabel("Cell count")
    ax_gen.set_title(f"Generations at t={END_TIME:.0f}h")
    ax_gen.set_xticks(unique_gens)
    ax_gen.grid(True, axis="y", color="#21262d", linewidth=0.5)

    # ── Animation ──
    all_pos  = np.vstack([pos for _, pos, _, _ in snapshots])
    pad      = 2.0
    xlim     = (all_pos[:, 0].min() - pad, all_pos[:, 0].max() + pad)
    ylim     = (all_pos[:, 1].min() - pad, all_pos[:, 1].max() + pad)

    ax_main.set_xlim(xlim)
    ax_main.set_ylim(ylim)
    ax_main.set_aspect("equal")
    ax_main.set_xlabel("x (cell diameters)")
    ax_main.set_ylabel("y (cell diameters)")
    ax_main.set_title("Off-Lattice Colony (Chaste-style)")

    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax_main, fraction=0.035, pad=0.02)
    cbar.set_label("Generation", color="#8b949e")
    cbar.ax.yaxis.set_tick_params(color="#8b949e")
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="#8b949e")

    time_text   = ax_main.text(0.02, 0.96, "", transform=ax_main.transAxes,
                               color="#e6edf3", fontsize=9, va="top")
    count_text  = ax_main.text(0.02, 0.91, "", transform=ax_main.transAxes,
                               color="#e6edf3", fontsize=9, va="top")
    patches: List[Circle] = []

    # vertical line tracking animation time on count plot
    vline = ax_count.axvline(x=0, color="#f78166", linewidth=1.2, alpha=0.8)

    def init():
        for p in patches:
            p.remove()
        patches.clear()
        return []

    def update(frame_idx):
        t, positions, generations, ages = snapshots[frame_idx]
        for p in patches:
            p.remove()
        patches.clear()

        for pos, gen, age in zip(positions, generations, ages):
            # Colour by generation; alpha encodes relative cell-cycle phase
            color = cmap(norm(gen))
            # Slightly shrink if very young (just divided)
            cycle_dur = MEAN_CELL_CYCLE_TIME
            growth    = min(1.0, age / (cycle_dur * 0.3))
            r         = CELL_RADIUS * (0.6 + 0.4 * growth)
            circle    = Circle(pos, radius=r, color=color,
                               alpha=0.85, linewidth=0.4,
                               edgecolor="white")
            ax_main.add_patch(circle)
            patches.append(circle)

        time_text.set_text(f"t = {t:.1f} h")
        count_text.set_text(f"N = {len(positions)}")
        vline.set_xdata([t, t])
        return patches + [time_text, count_text, vline]

    ani = animation.FuncAnimation(
        fig, update, frames=len(snapshots),
        init_func=init, interval=60, blit=False
    )

    fig.suptitle("Chaste-Style Cell Growth & Division  |  Single-Cell Start",
                 color="#e6edf3", fontsize=13, fontweight="bold")

    return fig, ani


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("Running Chaste-style off-lattice simulation…")
    print(f"  End time : {END_TIME} h | dt = {DT} h | seed = {SEED}")

    snapshots, cell_count_ts = run_simulation()

    final_count = cell_count_ts[-1][1]
    print(f"  Final cell count : {final_count}")
    print(f"  Snapshots recorded : {len(snapshots)}")

    fig, ani = visualise(snapshots, cell_count_ts)

    out_gif = "/mnt/user-data/outputs/chaste_simulation.gif"
    out_png = "/mnt/user-data/outputs/chaste_final_frame.png"

    print("Saving animation (GIF) …")
    ani.save(out_gif, writer="pillow", fps=18, dpi=110)
    print(f"  → {out_gif}")

    # Save final frame as static PNG too
    # Advance to last frame manually
    t, positions, generations, ages = snapshots[-1]
    ax_main = fig.axes[0]
    for patch in list(ax_main.patches):
        patch.remove()
    import matplotlib.cm as mcm
    cmap2 = mcm.plasma
    norm2 = Normalize(vmin=0, vmax=max(generations.max(), 1))
    for pos, gen, age in zip(positions, generations, ages):
        color  = cmap2(norm2(gen))
        circle = plt.Circle(pos, radius=CELL_RADIUS, color=color,
                             alpha=0.85, linewidth=0.4, edgecolor="white")
        ax_main.add_patch(circle)

    fig.savefig(out_png, dpi=140, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    print(f"  → {out_png}")
    print("Done.")
