"""
Chaste-style interactive cell growth & division simulation.
Single window: colony view only, with Play/Pause, Reset, Speed slider.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.widgets import Button, Slider
from matplotlib.animation import FuncAnimation

# ── Parameters ────────────────────────────────────────────────────────────────
CELL_RADIUS    = 0.5
SPRING_K       = 15.0
DAMPING        = 1.0
CUTOFF         = 1.5 * 2 * CELL_RADIUS
DT             = 0.005
MEAN_CYCLE     = 4.0    # hours
SD_CYCLE       = 0.5
BIRTH_SEP      = 0.1
MAX_CELLS      = 80
SEED           = 42

np.random.seed(SEED)

GEN_CMAP = plt.cm.plasma

# ── Cell data (pure numpy/lists, no dataclasses) ───────────────────────────
_id = 0

def new_cell(x, y, generation, birth_time):
    global _id
    c = dict(
        id=_id, x=x, y=y,
        gen=generation,
        birth=birth_time,
        duration=max(2.0, np.random.normal(MEAN_CYCLE, SD_CYCLE)),
    )
    _id += 1
    return c

def compute_forces(cells):
    n = len(cells)
    fx = np.zeros(n); fy = np.zeros(n)
    xs = np.array([c['x'] for c in cells])
    ys = np.array([c['y'] for c in cells])
    for i in range(n):
        dx = xs - xs[i]; dy = ys - ys[i]
        dist = np.sqrt(dx*dx + dy*dy)
        dist[i] = np.inf
        mask = dist < CUTOFF
        if not mask.any():
            continue
        nat = 2 * CELL_RADIUS
        f = np.where(mask, SPRING_K * (dist - nat) / np.maximum(dist, 1e-9), 0.0)
        fx[i] = (f * dx)[mask].sum()
        fy[i] = (f * dy)[mask].sum()
    return fx, fy

# ── Simulation state ──────────────────────────────────────────────────────────
class Sim:
    def __init__(self):
        self.reset()

    def reset(self):
        global _id
        _id = 0
        np.random.seed(SEED)
        self.t = 0.0
        self.cells = [new_cell(0.0, 0.0, 0, 0.0)]
        self.done = False

    def step(self, n_steps=1):
        if self.done:
            return
        for _ in range(n_steps):
            self.t += DT
            # Divide
            to_div = [c for c in self.cells if (self.t - c['birth']) >= c['duration']]
            if to_div and len(self.cells) < MAX_CELLS:
                div_ids = {c['id'] for c in to_div}
                new = []
                for p in to_div:
                    if len(self.cells) - len(div_ids) + len(new) >= MAX_CELLS:
                        break
                    a = np.random.uniform(0, 2*np.pi)
                    dx, dy = np.cos(a)*BIRTH_SEP, np.sin(a)*BIRTH_SEP
                    new.append(new_cell(p['x']+dx, p['y']+dy, p['gen']+1, self.t))
                    new.append(new_cell(p['x']-dx, p['y']-dy, p['gen']+1, self.t))
                self.cells = [c for c in self.cells if c['id'] not in div_ids] + new
            # Forces + integrate
            fx, fy = compute_forces(self.cells)
            for i, c in enumerate(self.cells):
                c['x'] += DT / DAMPING * fx[i]
                c['y'] += DT / DAMPING * fy[i]
            if len(self.cells) >= MAX_CELLS:
                self.done = True
                break

# ── Build figure ──────────────────────────────────────────────────────────────
sim = Sim()
running = [False]

fig = plt.figure(figsize=(7, 8), facecolor='#0d1117')
fig.canvas.manager.set_window_title('Chaste Cell Division')

# Main axes (leaves room at bottom for controls)
ax = fig.add_axes([0.05, 0.18, 0.90, 0.78])
ax.set_facecolor('#0d1117')
ax.set_aspect('equal')
ax.tick_params(colors='#8b949e')
for spine in ax.spines.values():
    spine.set_edgecolor('#30363d')
ax.set_xlabel('x (cell diameters)', color='#8b949e', fontsize=9)
ax.set_ylabel('y (cell diameters)', color='#8b949e', fontsize=9)
ax.set_title('Chaste · Off-Lattice Colony', color='#e6edf3', fontsize=11, pad=8)

VPAD = 14
ax.set_xlim(-VPAD, VPAD)
ax.set_ylim(-VPAD, VPAD)

# Scatter plot handle
max_gen = 8
norm = plt.Normalize(vmin=0, vmax=max_gen)
scat = ax.scatter([], [], s=(2*CELL_RADIUS*72*fig.get_figwidth()/28)**2,
                  c=[], cmap=GEN_CMAP, norm=norm,
                  edgecolors='white', linewidths=0.3, alpha=0.9, zorder=3)

# Cycle-progress rings (one artist per cell, rebuilt each frame)
ring_artists = []

# Colorbar
sm = plt.cm.ScalarMappable(cmap=GEN_CMAP, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.01)
cbar.set_label('Generation', color='#8b949e', fontsize=8)
cbar.ax.tick_params(colors='#8b949e', labelsize=7)

info_text = ax.text(0.02, 0.97, '', transform=ax.transAxes,
                    color='#e6edf3', fontsize=8.5, va='top',
                    fontfamily='monospace',
                    bbox=dict(facecolor='#161b22', edgecolor='#30363d',
                              boxstyle='round,pad=0.4', alpha=0.8))

done_text = ax.text(0.5, 0.5, '', transform=ax.transAxes,
                    color='#ffd93d', fontsize=13, ha='center', va='center',
                    fontweight='bold',
                    bbox=dict(facecolor='#0d1117', edgecolor='#ffd93d',
                              boxstyle='round,pad=0.6', alpha=0.9))

# ── Widget axes ───────────────────────────────────────────────────────────────
ax_play  = fig.add_axes([0.05, 0.09, 0.18, 0.055])
ax_reset = fig.add_axes([0.27, 0.09, 0.18, 0.055])
ax_speed = fig.add_axes([0.54, 0.10, 0.40, 0.030])

btn_play  = Button(ax_play,  '▶  Play',  color='#161b22', hovercolor='#21262d')
btn_reset = Button(ax_reset, '↺  Reset', color='#161b22', hovercolor='#21262d')
slider    = Slider(ax_speed, 'Speed', 1, 20, valinit=4, valstep=1,
                   color='#58a6ff', initcolor='#58a6ff')

for btn in (btn_play, btn_reset):
    btn.label.set_color('#e6edf3')
    btn.label.set_fontsize(9)
slider.label.set_color('#8b949e')
slider.label.set_fontsize(8)
slider.valtext.set_color('#8b949e')
ax_speed.set_facecolor('#161b22')

# ── Draw frame ────────────────────────────────────────────────────────────────
def redraw():
    cells = sim.cells
    if not cells:
        return

    xs  = np.array([c['x'] for c in cells])
    ys  = np.array([c['y'] for c in cells])
    gens = np.array([c['gen'] for c in cells])
    ages = np.array([sim.t - c['birth'] for c in cells])
    durs = np.array([c['duration'] for c in cells])
    phases = np.clip(ages / durs, 0, 1)

    scat.set_offsets(np.c_[xs, ys])
    scat.set_array(gens.astype(float))

    # Remove old rings
    for a in ring_artists:
        a.remove()
    ring_artists.clear()

    # Draw cycle-progress arcs
    for i, c in enumerate(cells):
        ph = phases[i]
        color = GEN_CMAP(norm(c['gen']))
        lw   = 2.0 if ph > 0.85 else 1.0
        ec   = 'white' if ph > 0.85 else (*color[:3], 0.55)
        arc = mpatches.Arc(
            (c['x'], c['y']),
            width=2*CELL_RADIUS*1.45, height=2*CELL_RADIUS*1.45,
            angle=-90, theta1=0, theta2=ph*360,
            color=ec, linewidth=lw, zorder=4
        )
        ax.add_patch(arc)
        ring_artists.append(arc)

    info_text.set_text(
        f"t = {sim.t:.2f} h\n"
        f"cells = {len(cells)}\n"
        f"max gen = {max(gens)}"
    )

    if sim.done:
        done_text.set_text(f'Max cells reached ({MAX_CELLS})')
    else:
        done_text.set_text('')

    fig.canvas.draw_idle()

# ── Animation ─────────────────────────────────────────────────────────────────
def animate(_frame):
    if running[0] and not sim.done:
        n = int(slider.val)
        sim.step(n)
        redraw()

ani = FuncAnimation(fig, animate, interval=30, cache_frame_data=False)

# ── Button callbacks ──────────────────────────────────────────────────────────
def toggle_play(event):
    if sim.done:
        return
    running[0] = not running[0]
    btn_play.label.set_text('⏸  Pause' if running[0] else '▶  Play')
    fig.canvas.draw_idle()

def do_reset(event):
    running[0] = False
    btn_play.label.set_text('▶  Play')
    sim.reset()
    redraw()

btn_play.on_clicked(toggle_play)
btn_reset.on_clicked(do_reset)

# Initial draw
redraw()
plt.savefig('/mnt/user-data/outputs/chaste_interactive_preview.png',
            dpi=130, bbox_inches='tight', facecolor=fig.get_facecolor())
plt.show()
