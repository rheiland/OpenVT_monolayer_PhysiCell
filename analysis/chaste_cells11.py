"""
Reproduces Chaste's RepulsionForce / GeneralisedLinearSpringForce
CalculateForceBetweenNodes logic for a NodeBasedCellPopulation.

Force on node A due to node B (vector points from A toward B):

  overlap = distance - rest_length        (rest_length = r_a + r_b)

  if overlap < 0  (cells overlap -> repulsion):
      F_A = mu * log(1 + overlap / rest_length) * unit_AB
      (log argument is in (0,1), so log < 0, force points away from B)

  if overlap >= 0 (cells separated -> weak attraction, decays exponentially):
      F_A = mu * overlap * exp(-alpha * overlap / rest_length) * unit_AB
      (RepulsionForce never reaches this branch because it only calls
       CalculateForceBetweenNodes when distance < rest_length)

Parameters matching Relax11.cpp:
  mu    = 5.0   (SetMeinekeSpringStiffness)
  alpha = 5.0   (hard-coded in GeneralisedLinearSpringForce.cpp)
"""

import numpy as np


def calculate_force_between_nodes(pos_a, pos_b, r_a, r_b, mu=5.0, alpha=5.0):
    """
    Return force vector on node A due to node B.

    Parameters
    ----------
    pos_a, pos_b : array-like, shape (2,)
    r_a, r_b     : float  cell radii
    mu           : float  spring stiffness (MeinekeSpringStiffness)
    alpha        : float  exponential decay constant (fixed at 5.0 in Chaste)

    Returns
    -------
    force_a : ndarray shape (2,)  force on A  (force on B = -force_a)
    """
    pos_a = np.asarray(pos_a, dtype=float)
    pos_b = np.asarray(pos_b, dtype=float)

    diff = pos_b - pos_a          # vector from A to B
    distance = np.linalg.norm(diff)
    assert distance > 0, "Nodes must be at distinct positions"

    unit_ab = diff / distance
    rest_length = r_a + r_b
    overlap = distance - rest_length

    if overlap < 0:
        # Repulsion: logarithmic (matches Chaste's "stable" form)
        magnitude = mu * np.log(1.0 + overlap / rest_length)
    else:
        # Attraction: linear * exponential decay
        magnitude = mu * overlap * np.exp(-alpha * overlap / rest_length)

    return magnitude * unit_ab


def repulsion_force_contribution(positions, radii, mu=5.0, alpha=5.0):
    """
    Compute net force on every node, summing over all overlapping pairs
    (mirrors RepulsionForce::AddForceContribution).

    Parameters
    ----------
    positions : array-like, shape (N, 2)
    radii     : array-like, shape (N,)
    mu        : float
    alpha     : float

    Returns
    -------
    forces : ndarray shape (N, 2)
    """
    positions = np.asarray(positions, dtype=float)
    radii = np.asarray(radii, dtype=float)
    n = len(positions)
    forces = np.zeros_like(positions)

    for i in range(n):
        for j in range(i + 1, n):
            rest_length = radii[i] + radii[j]
            diff = positions[j] - positions[i]
            distance = np.linalg.norm(diff)
            if distance < rest_length:          # only repulsion, matching RepulsionForce
                f_i = calculate_force_between_nodes(
                    positions[i], positions[j], radii[i], radii[j], mu, alpha
                )
                forces[i] += f_i
                forces[j] -= f_i

    return forces


# ---------------------------------------------------------------------------
# Forward-Euler simulation loop matching Relax11.cpp
# ---------------------------------------------------------------------------
def simulate(positions, radii, mu=5.0, dt=0.01, end_time=5.0,
             sampling_multiple=10, damping=1.0):
    """
    Evolve positions using Forward Euler, matching Chaste's
    ForwardEulerNumericalMethod:
        new_pos = old_pos + dt * force / damping

    Chaste's default damping constant (mDampingConstantNormal) is 1.0.

    Parameters
    ----------
    positions        : ndarray (N, 2)  initial positions
    radii            : ndarray (N,)
    mu               : float  spring stiffness
    dt               : float  timestep
    end_time         : float
    sampling_multiple: int    record every this many steps (matches
                              SetSamplingTimestepMultiple in Relax11.cpp)
    damping          : float  damping constant

    Returns
    -------
    records : list of (time, positions_copy) tuples at sampled steps
    """
    positions = positions.copy()
    n_steps   = round(end_time / dt)
    records   = []

    for step in range(n_steps + 1):
        t = step * dt
        if step % sampling_multiple == 0:
            records.append((t, positions.copy()))

        if step < n_steps:
            forces    = repulsion_force_contribution(positions, radii, mu=mu)
            positions = positions + dt * forces / damping

    return records


# ---------------------------------------------------------------------------
# Reproduce the Relax11 initial configuration and run the full simulation
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import csv, sys

    N        = 11
    RADIUS   = 5.0
    OVERLAP  = 5.0
    SPACING  = 2 * RADIUS - OVERLAP   # = 5.0
    MU       = 5.0
    DT       = 0.01
    END_TIME = 5.0
    SAMPLING = 10
    CSV_OUT  = "cells_py.csv"

    positions = np.array([[(i - 5) * SPACING, 0.0] for i in range(N)])
    radii     = np.full(N, RADIUS)

    print("Initial positions (x only, y=0):")
    print(" ", positions[:, 0])

    records = simulate(positions, radii, mu=MU, dt=DT, end_time=END_TIME,
                       sampling_multiple=SAMPLING)

    # Write CSV matching CsvOutputModifier: time,cell_id,x,y,radius
    with open(CSV_OUT, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["time", "cell_id", "x", "y", "radius"])
        for t, pos in records:
            for i in range(N):
                writer.writerow([f"{t:.4f}", i,
                                  f"{pos[i,0]:.6f}", f"{pos[i,1]:.6f}",
                                  RADIUS])

    print(f"\nWrote {len(records)} sampled timesteps to {CSV_OUT}")

    # Print final positions
    t_final, pos_final = records[-1]
    print(f"\nFinal positions at t={t_final:.2f} (x only):")
    for i in range(N):
        print(f"  cell {i:2d}  x={pos_final[i,0]:+.6f}")

