import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt

# Example points
# points = np.array([[0, 0], [0, 2], [1, 1], [2, 0], [2, 2]])
#points = np.loadtxt('pts_378.csv', delimiter=',')
points = np.loadtxt('run01/cells.csv', delimiter=',')

# Compute the convex hull
hull = ConvexHull(points)

# Get the vertices of the convex hull
hull_vertices = points[hull.vertices]

# Plotting (optional)
plt.plot(points[:, 0], points[:, 1], 'o', label='Points') # Plot the original points
plt.plot(hull_vertices[:, 0], hull_vertices[:, 1], 'r-', label='Convex Hull') # Plot the hull
plt.fill(hull_vertices[:, 0], hull_vertices[:, 1], 'r', alpha=0.3) # Fill the hull
plt.legend()
# plt.show()

# print("Vertices of the convex hull:", hull_vertices)
print("# of Vertices of the convex hull:", len(hull_vertices))

convex_hull_ids = []
for idx in range(len(hull_vertices)):
    for jdx in range(len(points)):
        if hull_vertices[idx][0] == points[jdx][0]:
            print(f'{idx}: {hull_vertices[idx][0]} at {jdx}')
            convex_hull_ids.append(jdx)

print(convex_hull_ids)
