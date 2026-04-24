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
# Reproduce the Relax11 initial configuration and step forward one timestep
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    N       = 11
    RADIUS  = 5.0
    OVERLAP = 5.0
    SPACING = 2 * RADIUS - OVERLAP   # = 5.0
    MU      = 5.0
    DT      = 0.01

    positions = np.array([[(i - 5) * SPACING, 0.0] for i in range(N)])
    radii     = np.full(N, RADIUS)

    print("Initial positions (x only, y=0):")
    print("  ", positions[:, 0])

    forces = repulsion_force_contribution(positions, radii, mu=MU)

    print("\nForces after one evaluation (x component):")
    for i in range(N):
        print(f"  cell {i:2d}  fx={forces[i, 0]:+.6f}  fy={forces[i, 1]:+.6f}")

    # Forward-Euler step  (Chaste uses damped mechanics: velocity = force / damping)
    # Chaste's default damping constant is 1.0, so displacement = force * dt
    new_positions = positions + forces * DT
    print("\nPositions after one Euler step (x only):")
    print("  ", new_positions[:, 0])
