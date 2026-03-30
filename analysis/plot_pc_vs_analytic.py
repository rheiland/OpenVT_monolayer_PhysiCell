import numpy as np
import matplotlib.pyplot as plt

#x, y = np.loadtxt('pc_plot_11cells.csv', delimiter=',', unpack=True)
x, y = np.loadtxt('pc_plot_11cells.csv', delimiter=',', skiprows=1, unpack=True)
plt.plot(x, y,'k',label='PhysiCell')
plt.legend()

x, y = np.loadtxt('t_dist_11cells_analytic.csv', delimiter=',',  unpack=True)
plt.plot(x, y,color='blue',linestyle='--',label='analytic')
plt.legend()

plt.plot([1,1], [5,10],'g--')
plt.plot([0,10], [9,9],'g--')

plt.xlabel('time')
plt.ylabel('distance between outer cells')
plt.title('11 cells relaxation')
plt.show()

