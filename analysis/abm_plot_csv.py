"""
plot_abm.py
-----------
Visualise agent-based model (ABM) data from abm_AB.csv.

Each agent is drawn as a filled circle (radius = radius) with an arrow
showing its instantaneous velocity vector (derived from successive positions).
Circles are coloured by norm_rand_i.

Usage
-----
    python plot_abm.py                      # interactive slider
    python plot_abm.py --time 1000          # single snapshot saved to PNG
    python plot_abm.py --animate            # save MP4 animation
    python plot_abm.py --data /path/to.csv  # custom CSV path
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.widgets import Slider
import matplotlib.cm as cm
import matplotlib.colors as mcolors

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
DEFAULT_CSV = "abm_AB.csv"
CMAP_NAME   = "plasma"       # colormap for norm_rand_i
ALPHA_CIRC  = 0.65           # circle fill opacity
ALPHA_EDGE  = 0.9            # circle edge opacity
ARROW_SCALE = 5.0            # length multiplier for velocity arrows
ARROW_WIDTH = 0.3            # quiver arrow width
FIG_SIZE    = (9, 8)


# ---------------------------------------------------------------------------
# Data loading & preprocessing
# ---------------------------------------------------------------------------
def load_data(csv_path: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df = df.sort_values(["ID", "time"]).reset_index(drop=True)

    # Compute velocity from position differences (x_vec/y_vec are stored as 0)
    df["vx"] = df.groupby("ID")["x_pos"].diff().fillna(0.0)
    df["vy"] = df.groupby("ID")["y_pos"].diff().fillna(0.0)

    # If the CSV already has non-zero vectors, prefer those
    if (df["x_vec"] != 0).any():
        df["vx"] = df["x_vec"]
        df["vy"] = df["y_vec"]

    return df


def get_frame(df: pd.DataFrame, t: float) -> pd.DataFrame:
    """Return all agents alive at time t."""
    return df[df["time"] == t].copy()


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------
def draw_frame(ax, frame: pd.DataFrame, cmap, norm, show_ids: bool = True):
    """
    Draw one simulation frame onto *ax*.
    Returns the list of artists added so they can be cleared next frame.
    """
    ax.cla()

    if frame.empty:
        ax.set_title("No data for this timestep")
        return

    t_val = frame["time"].iloc[0]

    # --- circles ------------------------------------------------------------
    circles = []
    colors  = []
    for _, row in frame.iterrows():
        circ = mpatches.Circle((row.x_pos, row.y_pos), row.radius)
        circles.append(circ)
        colors.append(norm(row.norm_rand_i))

    pc = PatchCollection(
        circles,
        cmap=cmap,
        norm=mcolors.Normalize(vmin=0, vmax=1),
        alpha=ALPHA_CIRC,
        edgecolors="white",
        linewidths=0.8,
    )
    pc.set_array(np.array(colors))
    ax.add_collection(pc)

    # --- velocity arrows ----------------------------------------------------
    has_motion = (frame["vx"] != 0) | (frame["vy"] != 0)
    moving = frame[has_motion]
    if not moving.empty:
        ax.quiver(
            moving["x_pos"].values,
            moving["y_pos"].values,
            moving["vx"].values * ARROW_SCALE,
            moving["vy"].values * ARROW_SCALE,
            color="white",
            alpha=0.9,
            width=ARROW_WIDTH,
            scale=1,
            scale_units="xy",
            angles="xy",
            zorder=5,
        )

    # --- agent ID labels ----------------------------------------------------
    if show_ids and len(frame) <= 50:
        for _, row in frame.iterrows():
            ax.text(
                row.x_pos, row.y_pos,
                str(int(row.ID)),
                ha="center", va="center",
                fontsize=6, color="white", fontweight="bold", zorder=6,
            )

    # --- axes / title -------------------------------------------------------
    pad = 15
    all_x = df_global["x_pos"]
    all_y = df_global["y_pos"]
    ax.set_xlim(all_x.min() - pad, all_x.max() + pad)
    ax.set_ylim(all_y.min() - pad, all_y.max() + pad)
    ax.set_aspect("equal")
    ax.set_facecolor("#1a1a2e")
    ax.set_xlabel("x position", color="white")
    ax.set_ylabel("y position", color="white")
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_edgecolor("#444")
    ax.set_title(
        f"ABM Snapshot  |  t = {t_val:.0f}  |  agents = {len(frame)}",
        color="white", fontsize=12, pad=10,
    )


# ---------------------------------------------------------------------------
# Single snapshot
# ---------------------------------------------------------------------------
def plot_snapshot(df: pd.DataFrame, t: float, save_path: str = None):
    # cmap = cm.get_cmap(CMAP_NAME)
    cmap = mpatplotlib.colormaps.get_cmap(CMAP_NAME)
    vmin, vmax = df["norm_rand_i"].min(), df["norm_rand_i"].max()
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    global df_global
    df_global = df

    fig, ax = plt.subplots(figsize=FIG_SIZE, facecolor="#0f0f1e")
    frame = get_frame(df, t)
    draw_frame(ax, frame, cmap, norm)

    # colorbar
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("norm_rand_i", color="white")
    cbar.ax.yaxis.set_tick_params(color="white")
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="white")

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
        print(f"Saved → {save_path}")
    else:
        plt.show()
    plt.close(fig)


# ---------------------------------------------------------------------------
# Interactive slider
# ---------------------------------------------------------------------------
def interactive_slider(df: pd.DataFrame):
    times = sorted(df["time"].unique())
    cmap  = cm.get_cmap(CMAP_NAME)
    vmin, vmax = df["norm_rand_i"].min(), df["norm_rand_i"].max()
    norm  = mcolors.Normalize(vmin=vmin, vmax=vmax)

    global df_global
    df_global = df

    fig = plt.figure(figsize=(FIG_SIZE[0], FIG_SIZE[1] + 1), facecolor="#0f0f1e")
    ax  = fig.add_axes([0.08, 0.12, 0.82, 0.80], facecolor="#1a1a2e")

    # colorbar axis
    cax = fig.add_axes([0.92, 0.12, 0.02, 0.80])
    sm  = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label("norm_rand_i", color="white")
    cbar.ax.yaxis.set_tick_params(color="white")
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="white")

    # slider
    ax_slider = fig.add_axes([0.12, 0.03, 0.75, 0.03], facecolor="#333355")
    slider = Slider(
        ax_slider, "Time", times[0], times[-1],
        valinit=times[0], valstep=times,
        color="#6666cc",
    )
    slider.label.set_color("white")
    slider.valtext.set_color("white")

    def update(val):
        t = slider.val
        frame = get_frame(df, t)
        draw_frame(ax, frame, cmap, norm)
        fig.canvas.draw_idle()

    slider.on_changed(update)
    update(times[0])

    plt.show()


# ---------------------------------------------------------------------------
# Animation
# ---------------------------------------------------------------------------
def save_animation(df: pd.DataFrame, out_path: str, fps: int = 15, step: int = 10):
    from matplotlib.animation import FFMpegWriter, FuncAnimation

    times = sorted(df["time"].unique())[::step]
    cmap  = cm.get_cmap(CMAP_NAME)
    vmin, vmax = df["norm_rand_i"].min(), df["norm_rand_i"].max()
    norm  = mcolors.Normalize(vmin=vmin, vmax=vmax)

    global df_global
    df_global = df

    fig, ax = plt.subplots(figsize=FIG_SIZE, facecolor="#0f0f1e")
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("norm_rand_i", color="white")
    cbar.ax.yaxis.set_tick_params(color="white")
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="white")

    def animate(t):
        frame = get_frame(df, t)
        draw_frame(ax, frame, cmap, norm)

    ani = FuncAnimation(fig, animate, frames=times, interval=1000 / fps)
    writer = FFMpegWriter(fps=fps, metadata={"title": "ABM simulation"})
    ani.save(out_path, writer=writer, dpi=120)
    print(f"Animation saved → {out_path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot ABM circles + velocity arrows.")
    parser.add_argument("--data",    default=DEFAULT_CSV, help="Path to CSV file")
    parser.add_argument("--time",    type=float,          help="Single timestep to plot")
    parser.add_argument("--save",    default=None,        help="Output PNG path (for --time)")
    parser.add_argument("--animate", action="store_true", help="Export MP4 animation")
    parser.add_argument("--anim-out",default="abm_animation.mp4", help="Animation output path")
    parser.add_argument("--fps",     type=int, default=15,  help="Frames per second")
    parser.add_argument("--step",    type=int, default=10,  help="Timestep stride for animation")
    args = parser.parse_args()

    print(f"Loading data from: {args.data}")
    df = load_data(args.data)
    print(f"  {len(df):,} rows | {df['ID'].nunique()} unique IDs | "
          f"t = {df['time'].min():.0f}–{df['time'].max():.0f}")

    if args.animate:
        save_animation(df, args.anim_out, fps=args.fps, step=args.step)

    elif args.time is not None:
        available = sorted(df["time"].unique())
        if args.time not in available:
            # snap to nearest
            t = min(available, key=lambda x: abs(x - args.time))
            print(f"  t={args.time} not found, using nearest t={t}")
        else:
            t = args.time
        plot_snapshot(df, t, save_path=args.save)

    else:
        # Default: interactive slider
        print("Launching interactive slider (close window to exit)…")
        interactive_slider(df)
