import numpy as np
import matplotlib.pyplot as plt

# Initial conditions
ics = np.arange(0, 5.5, 0.5)

# Time array
times = np.linspace(0, 10, 100)

start_time = times[0]

num_cells = len(ics)

# Initialize equations matrix
eqns = np.zeros((num_cells, num_cells))

k = 18.2816648

# Set up equations matrix
eqns[0, -1] = 1
eqns[1, -2:] = [1/k, 1]

for row in range(2, num_cells):
    eqns[row, :] = np.concatenate(([0], 1/k * eqns[row-1, 1:])) + 2*eqns[row-1, :] - eqns[row-2, :]

# Construct polynomial
p = np.concatenate((1/k * eqns[-1, :], [0])) + np.concatenate(([0], eqns[-1, :])) - np.concatenate(([0], eqns[-2, :]))

# Find roots (numpy.roots expects highest degree first, opposite of MATLAB)
L = np.roots(p)

# Build matrices
L_pow = np.zeros((num_cells, len(L)))

for row in range(num_cells):
    L_pow[num_cells - 1 - row, :] = L ** row

A = np.zeros((num_cells, len(L)))
for row in range(num_cells):
    A[row, :] = eqns[row, :] @ L_pow

# Solve for coefficients
b = -(np.arange(1, num_cells + 1)) + ics

coefficients = np.linalg.solve(A, b)

# Calculate positions over time
r = np.zeros((num_cells, len(times)))

for row in range(num_cells):
    r[row, :] = (eqns[row, :] @ L_pow * coefficients) @ np.exp(np.outer(L, times - start_time)) + (row + 1)

# Calculate separation
sep = ((eqns[-1, :] @ L_pow * coefficients) @ np.exp(np.outer(L, times - start_time)) + 11) - \
      ((eqns[0, :] @ L_pow * coefficients) @ np.exp(np.outer(L, times - start_time)) + 1)

# Plot positions
plt.figure()
plt.plot(times, r.T)
plt.xlabel('Time')
plt.ylabel('Position')
plt.title('Spring Positions Over Time')
plt.grid(True)

# Plot separation
plt.figure()
plt.plot(times, sep, label='Separation')
plt.axhline(y=9, color='k', linestyle=':', label='Reference line (y=9)')
plt.axvline(x=1, color='k', linestyle=':', label='Reference line (x=1)')
plt.xlabel('Time')
plt.ylabel('Separation')
plt.title('Separation Between First and Last Spring')
plt.legend()
plt.grid(True)

# Calculate separation at t=1
sep_1 = ((eqns[-1, :] @ L_pow * coefficients) @ np.exp(L) + 11) - \
        ((eqns[0, :] @ L_pow * coefficients) @ np.exp(L) + 1)

print(f"Separation at t=1: {sep_1}")

plt.show()