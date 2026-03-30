"""
Growing monolayer model in Python 
=================================

This initial version has no contact inhibition. The goal is to understand why PhysiCell seems
to differ from Chaste, in terms of overcrowding/overlapping cells in the center of the monolayer
when it gets too large.

* Cell size maintenance (parameter: target cell area A_0(t))
* start with a single cell at the origin, with area = A_0(t)
* Cell growth (parameters: linear growth rate α)
* Cell division: cells divide when A(t) = X * A_0(0) 
  (parameter:  X ∼ N(2, 0.4^2), truncated such that X > 0, 
   redraw X otherwise). 
   Observable for an isolated cell: cell cycle duration T = A0(0)/α)


This is the slower implementation of the model, to hopefully make it more obvious how
cells grow and divide. (The faster version is abm.py, uses numpy, and is a bit more challenging to see
the details.)

"""
import sys
import math
import random
# import uuid
import csv
from dataclasses import dataclass, field
from typing import Optional
import matplotlib.pyplot as plt
import matplotlib.patches as patches
# import mpld3

print("# args=",len(sys.argv))
if len(sys.argv) < 5:
    print("Usage: <repulsion(10)> <max cells> <win size> <split type(0-2)>\n")
    exit()
idx=1
repulsion = float(sys.argv[idx])
idx+=1
MAX_AGENTS = int(sys.argv[idx])
idx+=1
win_size = float(sys.argv[idx])
idx+=1
daughter_split = int(sys.argv[idx])

max_ID = 0

@dataclass
class Agent:
    """
    An agent with a 2D position, area, and constant growth rate.

    Attributes:
        x, y          : Position 
        vel_x, vel_y  : Velocity 
        area        : Current area of the agent.
        growth_rate   : Rate at which the area increases per time step (units/step).
        division_area : area at which the agent divides.
        norm_rand     : X ∼ N(2, 0.4^2), truncated such that X > 0, redraw X otherwise
        ID            : unique integer identifier
    """
    x: float
    y: float
    vel_x: float
    vel_y: float
    area: float = 78.5  # pi * 5^2
    growth_rate: float = 0.1778   # microns^2(?)/min  0.05
    division_area: float = 157.08        # division area=7.07
    norm_rand: float = 2.0
    ID: int = 0
    time_step: int = 0

    def grow(self) -> None:
        """Increase area by growth_rate for one time step."""
            #   -- may need to scale by dt here (Dom)
        self.area += self.growth_rate
        self.time_step += 1
        # if self.time_step % 5 == 0:

    def should_divide(self) -> bool:
        """Return True if this agent has reached or exceeded its division area."""
        return self.area >= self.division_area

    def divide(self, separation: float = 0.5, noise: float = 0.1) -> tuple["Agent", "Agent"]:
        global max_ID
        """
        Divide this agent into two daughter agents.

        The daughters have 1/2 the area of the parent. One daughter remains in the parent's position;
        the other daughter is placed at a random position around the other, so they just touch.

        Returns:
            A tuple of two new Agent instances (daughter_1, daughter_2).
        """

        daughter_area = self.area / 2
        daughter_radius = math.sqrt(daughter_area/math.pi)

        max_ID += 1

        theta = random.random() * 6.283185307179    # 2 * math.pi
        xvec = math.cos(theta) * daughter_radius * 2
        yvec = math.sin(theta) * daughter_radius * 2

        nr1 = random.normalvariate(mu=2.0, sigma=0.4)   # mean=2, stddev=0.4
        while nr1 < 0:
             nr1 = random.normalvariate(mu=2.0, sigma=0.4)   # mean=2, stddev=0.4

        nr2 = random.normalvariate(mu=2.0, sigma=0.4)   # mean=2, stddev=0.4
        while nr2 < 0:
             nr2 = random.normalvariate(mu=2.0, sigma=0.4)   # mean=2, stddev=0.4

        if daughter_split == 0:     # PhysiCell: daughter = parent x,y
            d1 = Agent(x=self.x,      y=self.y,      vel_x=0, vel_y=0,
                    area=daughter_area, growth_rate=self.growth_rate,
                    division_area=78.54 * nr1, norm_rand=nr1, ID=self.ID)
            d2 = Agent(x=self.x + xvec, y=self.y + yvec, vel_x=0, vel_y=0,
                    area=daughter_area, growth_rate=self.growth_rate,
                    division_area=78.54 * nr2, norm_rand=nr2, ID=max_ID)
        elif daughter_split == 1:   # equi-split
            xvec /= 2
            yvec /= 2
            d1 = Agent(x=self.x + xvec,      y=self.y + yvec,      vel_x=0, vel_y=0,
                    area=daughter_area, growth_rate=self.growth_rate,
                    division_area=78.54 * nr1, norm_rand=nr1, ID=self.ID)
            d2 = Agent(x=self.x - xvec, y=self.y - yvec, vel_x=0, vel_y=0,
                    area=daughter_area, growth_rate=self.growth_rate,
                    division_area=78.54 * nr2, norm_rand=nr2, ID=max_ID)
        elif daughter_split == 2:   # equi-split w/ overlap
            xvec /= 2.5
            yvec /= 2.5
            d1 = Agent(x=self.x + xvec,      y=self.y + yvec,      vel_x=0, vel_y=0,
                    area=daughter_area, growth_rate=self.growth_rate,
                    division_area=78.54 * nr1, norm_rand=nr1, ID=self.ID)
            d2 = Agent(x=self.x - xvec, y=self.y - yvec, vel_x=0, vel_y=0,
                    area=daughter_area, growth_rate=self.growth_rate,
                    division_area=78.54 * nr2, norm_rand=nr2, ID=max_ID)

        return d1, d2


class Simulation:
    """
    Manages a population of Agents over discrete time steps.

    At each step:
      1. Each cell grows
      2. Cells that exceed division_area are replaced by two daughters.
    """

    # Relevant parameter values will be passed in
    def __init__(
        self,
        n_seed: int = 0,
        growth_rate: float = 0.0,
        division_area: float = 0.0,
        norm_rand: float = 0.0,
        initial_area: float = 0.0,
    ):
        self.agents: list[Agent] = [
            Agent(
                x=0,   # position
                y=0,
                vel_x=0,   # velocity
                vel_y=0,
                area=initial_area,
                growth_rate=growth_rate,
                division_area=division_area,
                norm_rand=norm_rand,
            )
            for _ in range(n_seed)
        ]
        self.time: int = 0
        # self.history: list[dict] = []

    # involves the crucial assumptions related to "time"
    def step(self) -> int:
        """
        Advance the simulation by one time step.

        Returns:
            Number of division events that occurred this step.
        """
        survivors: list[Agent] = []
        divisions = 0

        # 1) Grow all agents
        for agent in self.agents:
            agent.grow()      # <----------------------
            if agent.should_divide():
                d1, d2 = agent.divide()
                survivors.extend([d1, d2])
                divisions += 1
            else:
                survivors.append(agent)
            # print("survivors=",survivors)

        self.agents = survivors

        # void standard_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt)
        # for(neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	    # {
		#     pCell->add_potentials(*neighbor);
	    # }

        # 2) Update each agent's velocity vector:  like Cell::add_potentials(Cell* other_agent)
        # cell_cell_repulsion_strength = 10.0
        cell_cell_repulsion_strength = repulsion
        # cell_cell_repulsion_strength = 0.0

        max_relax_steps = 1
        for relax_steps in range(max_relax_steps):   # isn't this violating our assumptions of time??

            for agent in self.agents:
                agent_r = math.sqrt(agent.area/math.pi)
                agent.vel_x = 0
                agent.vel_y = 0
                for agent2 in self.agents:
                    if agent.ID != agent2.ID:
                        # displacement[i] = position[i] - (*other_agent).position[i]; 
                        xdel = agent.x - agent2.x   # "displacement" vector in C++
                        ydel = agent.y - agent2.y
                        dist = math.sqrt(xdel*xdel + ydel*ydel)
                        # distance = std::max(sqrt(distance), 0.00001); 

        # //Repulsive
        # double R = phenotype.geometry.radius+ (*other_agent).phenotype.geometry.radius; 

                        agent2_r = math.sqrt(agent2.area/math.pi)
                        if dist < 1.25 * (agent_r + agent2_r):
                            R = agent_r + agent2_r
                            if dist > R:
                                temp_r = 0
                            else:
                                # print(f"update velocity of cell {agent.ID} due to nbr {agent.ID}")
                                temp_r = -dist  # -d
                                temp_r /= R # -d/R
                                temp_r += 1.0 # 1-d/R
                                temp_r *= temp_r # (1-d/R)^2 
                                
                                # add the relative pressure contribution 
                                # state.simple_pressure += ( temp_r / simple_pressure_scale ); // New July 2017 

                            # PhysiCell C++ code
                            # double effective_repulsion = sqrt( phenotype.mechanics.cell_cell_repulsion_strength * other_agent->phenotype.mechanics.cell_cell_repulsion_strength ); 
                            # effective_repulsion = math.sqrt(cell_cell_repulsion_strength * cell_cell_repulsion_strength) 
                            # effective_repulsion = cell_cell_repulsion_strength  
                            # temp_r *= effective_repulsion 
                            temp_r *= cell_cell_repulsion_strength 

                            # skip over adhesion, there is none

                            # if( fabs(temp_r) < 1e-16 )
                            # { return; }
                            temp_r /= dist

                            # velocity[i] += displacement[i] * temp_r; 
                            # axpy( &velocity , temp_r , displacement ); 
                            agent.vel_x += xdel * temp_r
                            agent.vel_y += ydel * temp_r
                            # if self.time >= 90 and self.time <=120:  #rwh
                                # print(f"--ID={agent.ID}, step={self.time}: update velocity for {agent.agent_id}: vel_x={agent.vel_x}, _y={agent.vel_y}")
                                # print(f"--step={self.time}: update velocity for {agent.ID}: vel_x={agent.vel_x}, _y={agent.vel_y}")

            # 3) Update each agent's position with its velocity
            #   -- may need to scale by dt here (Dom)
            for agent in self.agents:
                agent.x += agent.vel_x
                agent.y += agent.vel_y

        self.time += 1

        return divisions

    def run(self, steps: int, max_agents: int = 500) -> None:
        """
        Run the simulation for a given number of steps or until max_agents is reached.

        Args:
            steps:      Maximum number of time steps to run.
            max_agents: Stop early if population exceeds this threshold.
        """

        for _ in range(steps):
            self.step()
            if len(self.agents) >= max_agents:
                print(f"  Stopped early at t={self.time}: {len(self.agents)} agents.")
                break


    def precompute(self, steps: int, max_agents: int = 500) -> None:
        """
        Run the full simulation up-front and store every frame's agent list.
        Required before calling interactive_viewer().

        Args:
            steps:      Number of time steps to pre-compute.
            max_agents: Stop early if population exceeds this.
        """
        import copy
        print(f"Pre-computing {steps} steps...", end=" ", flush=True)
        # Frame 0 = initial state (before any stepping)
        self._frames: list[list[Agent]] = [copy.deepcopy(self.agents)]

        for _ in range(steps):
            self.step()
            self._frames.append(copy.deepcopy(self.agents))
            if len(self.agents) >= max_agents:
                print(f"(stopped early at t={self.time} — {len(self.agents)} agents)", end=" ")
                break

        print(f"done. {len(self._frames)} frames stored.")

        # save final results in .csv for post-processing, e.g., to compute f_i and a_i
        file_out = f'py_monolayer_small.csv'
        print("--> ",file_out)
        with open(file_out, "w", newline="") as file:
            writer = csv.writer(file)

            writer.writerow(['x_pos','y_pos','radius_i','norm_rand_i'])
            for agent in self.agents:
                # writer.writerow([x_pos[jdx],y_pos[jdx],radius_i[jdx],f_i[jdx],a_i[jdx]])
                radius = math.sqrt(agent.area / math.pi)
                writer.writerow([agent.x, agent.y, radius, agent.norm_rand])


    def interactive_viewer(self, title: str = "Agent-Based Model", interval: int = 100) -> None:
        """
        Open an interactive matplotlib window with:
          ◀◀  Step back one frame
          ▶▶  Step forward one frame
          ▶ / ■  Play / Pause the full simulation
          A time-step slider for direct scrubbing

        Requires precompute() to have been called first.
        """
        import copy
        from matplotlib.widgets import Button, Slider
        from matplotlib.animation import FuncAnimation

        if not hasattr(self, "_frames"):
            raise RuntimeError("Call precompute() before interactive_viewer().")

        frames     = self._frames
        n_frames   = len(frames)
        cmap       = plt.cm.plasma

        fig = plt.figure(figsize=(6, 6), facecolor="#0d0d0d")
        fig.suptitle(title, fontsize=13, fontweight="bold", color="white")

        ax_space = fig.add_axes([0.1, 0.2, 0.8, 0.7])   # left: spatial: left, bottom, w, h)

        ax_space.set_facecolor("#111111")
        ax_space.tick_params(colors="white")
        for sp in ax_space.spines.values():
            sp.set_edgecolor("#444444")

        ax_space.set_aspect("equal")
        ax_space.set_xlabel("x", color="white")
        ax_space.set_ylabel("y", color="white")
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 10))
        sm.set_array([])

        # Stable world bounds (union of all frame extents)
        all_x = [a.x for f in frames for a in f]
        all_y = [a.y for f in frames for a in f]
        # max_r = max((a.area for f in frames for a in f), default=1)  # now doing area, not radii
        # pad   = max_r * 2
        pad   = 20
        # ax_space.set_xlim(min(all_x) - pad, max(all_x) + pad)
        # ax_space.set_ylim(min(all_y) - pad, max(all_y) + pad)

        # win_size = 100  # rwh: hard-code for now
        # win_size = 150  # rwh: hard-code for now
        ax_space.set_xlim(-win_size, win_size)
        ax_space.set_ylim(-win_size, win_size)

        title_text = ax_space.set_title("", color="white", fontsize=10)

        # ------------------------------------------------------------------ #
        # Widget axes
        # ------------------------------------------------------------------ #
        btn_color      = "#1e2d3d"
        btn_hover      = "#2e4d6d"

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

        # ------------------------------------------------------------------ #
        # Draw helpers
        # ------------------------------------------------------------------ #
        state = {"frame": 0, "playing": False}

        def _draw_frame(idx: int) -> None:
            idx = max(0, min(idx, n_frames - 1))
            state["frame"] = idx

            agent_list = frames[idx]
            # max_gen    = max((a.generation for a in agent_list), default=0) or 1

            for p in list(ax_space.patches):
                p.remove()
            for agent in agent_list:
                # color = cmap(agent.generation / max(max_gen, 10))
                ax_space.add_patch(patches.Circle(
                    (agent.x, agent.y), 
                    radius=math.sqrt(agent.area / math.pi),
                    facecolor='gray', alpha=0.72,
                    linewidth=1.0, edgecolor="white",
                ))

            title_text.set_text(
                f"t = {idx}   |   agents = {len(agent_list)}"  #   |   max gen = {max_gen}"
            )

            # pop_line.set_data(all_times[:idx + 1], all_pops[:idx + 1])
            # time_marker.set_xdata([idx, idx])

            # Update slider without triggering its callback
            slider.eventson = False
            slider.set_val(idx)
            slider.eventson = True

            fig.canvas.draw_idle()

            # mpld3.save_html(fig, 'index.html') # Save directly to an HTML file



        # ------------------------------------------------------------------ #
        # Button / slider callbacks
        # ------------------------------------------------------------------ #
        def on_back(_event):
            state["playing"] = False
            btn_play.label.set_text("▶ Play")
            _draw_frame(state["frame"] - 1)

        def on_fwd(_event):
            state["playing"] = False
            btn_play.label.set_text("▶ Play")
            _draw_frame(state["frame"] + 1)

        def on_play(_event):
            state["playing"] = not state["playing"]
            btn_play.label.set_text("■ Pause" if state["playing"] else "▶ Play")
            fig.canvas.draw_idle()

        def on_slider(val):
            state["playing"] = False
            btn_play.label.set_text("▶ Play")
            _draw_frame(int(val))

        btn_back.on_clicked(on_back)
        btn_fwd.on_clicked(on_fwd)
        btn_play.on_clicked(on_play)
        slider.on_changed(on_slider)

        # ------------------------------------------------------------------ #
        # Timer-driven playback
        # ------------------------------------------------------------------ #
        def _tick(_frame):
            if state["playing"]:
                nxt = state["frame"] + 1
                if nxt >= n_frames:
                    state["playing"] = False
                    btn_play.label.set_text("▶ Play")
                else:
                    _draw_frame(nxt)

        anim = FuncAnimation(fig, _tick, interval=interval, cache_frame_data=False)

        # Draw initial frame
        _draw_frame(0)
        plt.show()

        return anim  # keep reference alive so GC doesn't kill the timer


# ---------------------------------------------------------------------------
# Example usage
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    from matplotlib.animation import FuncAnimation  # ensure import available

    random.seed(42)

    sim = Simulation(
        n_seed= 1,              # Start with n agents
        growth_rate= 0.1778,    # 0.889 original;  5x slower=0.1778  (micron^2 / min)
        division_area= 157.08,  #  vs. double area = 2 * 78.54 = 157.08
        norm_rand= 2.0,         #  X ~ N(2, 0.4^2)
        initial_area= 78.54,    #  area = pi * r^2 = pi * 25 = 78.54
    )

    # 1. Pre-compute all frames (enables backward scrubbing)
    sim.precompute(steps=50000, max_agents=MAX_AGENTS)  # rwh

    # 2. Open the interactive widget viewer
    sim.interactive_viewer(
        title="Growing Monolayer",
        interval=1,   # rwh: ms between auto-play frames
    )
