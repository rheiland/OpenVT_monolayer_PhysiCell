import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from matplotlib.colors import ListedColormap, BoundaryNorm
import csv

# %matplotlib widget


# import the needed data for each framework

def import_all_simulation_data(framework_head):
    # framework = '1000_Data/' +  framework_head +'_MonolayerGrowth_1000_Data/'
    framework = 'PhysiCell_MonolayerGrowth_1000_Data/'   # rwh

    # Aggregate data from all 100 simulations
    all_x_vals = []
    all_y_vals = []
    all_r_vals = []
    all_surface_fractions = []
    all_area_fractions = []
    all_inhibited_cells = []

    if framework_head == 'TST' or framework_head == 'Artistoo' or framework_head == 'Morpheus':
        iter_vals = range(1, 101)
    else:
        iter_vals = range(0, 100)

    for index in iter_vals:
        if framework_head == 'Artistoo':
            file_path = os.path.join(framework, f'cell_data_no_inhibition_{index:03d}.csv')
        else:
            file_path = os.path.join(framework, f'cell_data_no_inhibition_{index}.csv')
        
        if framework_head == 'Morpheus':
            data = np.loadtxt(file_path, delimiter='\t', skiprows=1)

            x_vals = data[:, 2]
            y_vals = data[:, 3]
            r_vals = data[:, 4]
            surface_fractions = data[:, 5]
            area_fractions = data[:, 6]
            
            all_x_vals.extend(x_vals)
            all_y_vals.extend(y_vals)
            all_r_vals.extend(r_vals)
            all_surface_fractions.extend(surface_fractions)
            all_area_fractions.extend(area_fractions)

        else:
            print("--- file_path= ",file_path)
            data = np.loadtxt(file_path, delimiter=',', skiprows=1)
            
            x_vals = data[:, 0]
            y_vals = data[:, 1]
            r_vals = data[:, 2]
            surface_fractions = data[:, 3]
            area_fractions = data[:, 4]
            
            all_x_vals.extend(x_vals)
            all_y_vals.extend(y_vals)
            all_r_vals.extend(r_vals)
            all_surface_fractions.extend(surface_fractions)
            all_area_fractions.extend(area_fractions)

    if framework_head == 'TST':
        tmp_sf = all_surface_fractions
        tmp_af = all_area_fractions

        all_surface_fractions = tmp_af
        all_area_fractions = tmp_sf


    all_surface_fractions = np.array(all_surface_fractions)
    all_surface_fractions[all_surface_fractions < 0.0] = 0.0

    all_area_fractions = np.array(all_area_fractions)
    all_area_fractions[all_area_fractions < 0.0] = 0.0

    return all_x_vals, all_y_vals, all_r_vals, all_surface_fractions, all_area_fractions


# framework_head = 'Artistoo'
# framework_head = 'CHASTE'
# framework_head = 'PolyHoop'
# framework_head = 'TST'
# framework_head = 'Morpheus'
framework_head = 'PhysiCell'

all_x_vals, all_y_vals, all_r_vals, all_surface_fractions, all_area_fractions = import_all_simulation_data(framework_head)

# compute distances and distance bins

y_0 = np.mean(all_y_vals)
x_0 = np.mean(all_x_vals)

distances = np.sqrt((np.array(all_x_vals) - x_0)**2 + (np.array(all_y_vals) - y_0)**2)
max_distance = np.ceil(np.max(distances))

# distances = distances/max_distance  # normalize distances to [0, 1]

number_of_bins = 7

bins = np.linspace(0, 1.05*np.max(distances), number_of_bins + 1)

def make_stacked(fractions, distances, bins):
    stacked = [[] for _ in range(len(bins) - 1)]
    for dist, val in zip(distances, fractions):
        idx = np.searchsorted(bins, dist, side='right') - 1
        if 0 <= idx < len(stacked):
            stacked[idx].append(val)
    return stacked

stacked_sf = make_stacked(all_surface_fractions, distances, bins)
stacked_af = make_stacked(all_area_fractions, distances, bins)

labels = [f'{int(bins[i])}-{int(bins[i+1])}' for i in range(len(bins)-1)]

# Plot stacked histograms
fig, axs = plt.subplots(2, 2, figsize=(10, 5))

# Surface fraction stacked histogram
xlims1 = (-0.01, 1.01)
cmap = plt.cm.viridis_r
axs[0, 0].hist(stacked_sf, bins=np.arange(xlims1[0], xlims1[1], 0.01), stacked=True, label=labels, edgecolor='none', color=cmap(np.linspace(0, 1, len(labels))), density=False)
axs[0, 0].set_title('PDF of Surface Fractions')
axs[0, 0].set_xlabel('Surface Fraction (γ)')
axs[0, 0].set_ylabel('Number of Cells')
axs[0, 0].set_yscale('log')
# axs[0].axvline(x=p_gamma, color='red', linestyle='--', label='γ Th.')
axs[0, 0].legend(title='Distance', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
# set axis limits
axs[0, 0].set_xlim(xlims1)
# axs[0].set_ylim(1, None)


xlims2 = (-0.01, 1.01)
cmap = plt.cm.inferno_r
# Area fraction stacked histogram
axs[1, 0].hist(stacked_af, bins=np.arange(xlims2[0], xlims2[1], 0.01), stacked=True, label=labels, edgecolor='none', color=cmap(np.linspace(0, 1, len(labels))), density=False)
axs[1, 0].set_title('PDF of Area Fractions')
axs[1, 0].set_xlabel('Area Fraction (β)')
axs[1, 0].set_ylabel('Number of Cells')
# axs[1, 0].set_yscale('log')
# axs[1].axvline(x=p_beta, color='red', linestyle='--', label='β Th.')
axs[1, 0].legend(title='Distance', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
axs[1, 0].set_xlim(xlims2)


# Surface fraction stacked histogram
# xlims1 = (0.7, 1.001)
cmap = plt.cm.viridis_r
axs[0, 1].hist(stacked_sf, bins=np.arange(xlims1[0], xlims1[1], 0.01), stacked=True, label=labels, edgecolor='none', color=cmap(np.linspace(0, 1, len(labels))), density=False, cumulative=True)
axs[0, 1].set_title('CDF of Surface Fractions')
axs[0, 1].set_xlabel('Surface Fraction (γ)')
axs[0, 1].set_ylabel('Number of Cells')
# axs[0].set_yscale('log')
# axs[0].axvline(x=p_gamma, color='red', linestyle='--', label='γ Th.')
axs[0, 1].legend(title='Distance', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
# set axis limits
axs[0, 1].set_xlim(xlims1)
# axs[0].set_ylim(1, None)


# xlims2 = (0.0, 1.01)
cmap = plt.cm.inferno_r
# Area fraction stacked histogram
axs[1, 1].hist(stacked_af, bins=np.arange(xlims2[0], xlims2[1], 0.01), stacked=True, label=labels, edgecolor='none', color=cmap(np.linspace(0, 1, len(labels))), density=False, cumulative=True)
axs[1, 1].set_title('CDF of Area Fractions')
axs[1, 1].set_xlabel('Area Fraction (β)')
axs[1, 1].set_ylabel('Number of Cells')
# axs[1].set_yscale('log')
# axs[1].axvline(x=p_beta, color='red', linestyle='--', label='β Th.')
axs[1, 1].legend(title='Distance', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
axs[1, 1].set_xlim(xlims2)

plt.suptitle(framework_head + ': Monolayer Growth - Surface and Area Fractions', fontweight='bold')
plt.tight_layout()
plt.show()


# fname = f'{framework_head}_MonolayerGrowth_Surface_Area_Fractions'
# fig.savefig('1000_Data/' + fname + '.png', dpi=300, bbox_inches='tight')

#------------------------------------------------

# framework_head = 'Artistoo'
# framework_head = 'CHASTE'
# framework_head = 'PolyHoop'
# framework_head = 'TST'
# framework_head = 'Morpheus'
framework_head = 'PhysiCell'
# framework = '1000_Data/' +  framework_head +'_MonolayerGrowth_1000_Data/'
framework = framework_head +'_MonolayerGrowth_1000_Data/'


# import the first simulation data for plotting
file_path = os.path.join(framework, f'cell_data_no_inhibition_1.csv')
data = np.loadtxt(file_path, delimiter=',', skiprows=1)
# data = np.loadtxt(file_path, delimiter='\t', skiprows=1)
x_vals = data[:, 0] - x_0
y_vals = data[:, 1] - y_0
r_vals = data[:, 2]
surface_fractions = data[:, 3]
surface_fractions[surface_fractions < 0.0] = 0.0
area_fractions = data[:, 4]
area_fractions[area_fractions < 0.0] = 0.0

p_beta = 0.0
p_gamma = 0.0

lower_val = -np.max(np.abs([x_vals, y_vals]))
upper_val = np.max(np.abs([x_vals, y_vals]))
Window_Lim_x = (lower_val, upper_val)
Window_Lim_y = (lower_val, upper_val)


# precompute normalizations for the colorbars
norm_r = plt.Normalize(min(r_vals), max(r_vals))
norm_sf = plt.Normalize(min(surface_fractions), 1)
norm_af = plt.Normalize(min(area_fractions), 1)
norm_ih = plt.Normalize(0, 3)

fig, axs = plt.subplots(1, 3, figsize=(12, 3))

# Plot 1: Colour by radius
axs[0].scatter(x_vals, y_vals, s=0.5, color='black')
for x, y, r in zip(x_vals, y_vals, r_vals):
    color = plt.cm.PRGn(norm_r(r))
    circle_1 = plt.Circle((x, y), r, color=color, fill=True, alpha=0.75)
    circle_2 = plt.Circle((x, y), r, color='black', fill=False, alpha=0.5)
    axs[0].add_artist(circle_1)
    axs[0].add_artist(circle_2)
axs[0].set_title('Cell Radius')
axs[0].set_xlabel('X Position')
axs[0].set_ylabel('Y Position')
axs[0].axis('equal')
axs[0].set_xlim(Window_Lim_x)
axs[0].set_ylim(Window_Lim_y)
# axs[0].grid(True)
sm_r = plt.cm.ScalarMappable(norm=norm_r, cmap='PRGn')
sm_r.set_array([])
fig.colorbar(sm_r, ax=axs[0], label='Cell Size (Radius)', ticks=np.linspace(min(r_vals), max(r_vals), 5))


#------------------------------------------------------------------------------#
# Plot 2: Colour by surface fraction
axs[1].scatter(x_vals, y_vals, s=0.5, color='black')
for x, y, sf, r in zip(x_vals, y_vals, surface_fractions, r_vals):
    color = plt.cm.viridis_r(norm_sf(sf))
    circle_1 = plt.Circle((x, y), r, color=color, fill=True, alpha=0.75)
    circle_2 = plt.Circle((x, y), r, color='black', fill=False, alpha=0.5)
    axs[1].add_artist(circle_1)
    axs[1].add_artist(circle_2)
axs[1].set_title('Surface Fraction, gamma')
axs[1].set_xlabel('X Position')
axs[1].set_ylabel('Y Position')
axs[1].axis('equal')
axs[1].set_xlim(Window_Lim_x)
axs[1].set_ylim(Window_Lim_y)
# axs[1].grid(True)
sm_sf = plt.cm.ScalarMappable(norm=norm_sf, cmap='viridis_r')
sm_sf.set_array([])
fig.colorbar(sm_sf, ax=axs[1], label='Surface Fraction')

#------------------------------------------------------------------------------#
# Plot 3: Colour by area fraction
axs[2].scatter(x_vals, y_vals, s=0.5, color='black')
for x, y, af, r in zip(x_vals, y_vals, area_fractions, r_vals):
    color = plt.cm.inferno_r(norm_af(af))
    circle_1 = plt.Circle((x, y), r, color=color, fill=True, alpha=0.75)
    circle_2 = plt.Circle((x, y), r, color='black', fill=False, alpha=0.5)
    axs[2].add_artist(circle_1)
    axs[2].add_artist(circle_2)
axs[2].set_title('Area Fraction, beta')
axs[2].set_xlabel('X Position')
axs[2].set_ylabel('Y Position')
axs[2].axis('equal')
axs[2].set_xlim(Window_Lim_x)
axs[2].set_ylim(Window_Lim_y)
# axs[2].grid(True)
sm_af = plt.cm.ScalarMappable(norm=norm_af, cmap='inferno_r')
sm_af.set_array([])
fig.colorbar(sm_af, ax=axs[2], label='Area Fraction')

plt.tight_layout()
plt.show()

#-----------------------------------------------------------------

plot_mean_std = True

# framework_heads = ['Artistoo', 'TST', 'Morpheus','PolyHoop', 'CHASTE', 'PhysiCell']
framework_heads = ['PhysiCell']

# Plot stacked histograms
fig, axs = plt.subplots(len(framework_heads), 4, figsize=(20, 13))


for iter, framework_head in enumerate(framework_heads):

    all_x_vals, all_y_vals, all_r_vals, all_surface_fractions, all_area_fractions = import_all_simulation_data(framework_head)


    y_0 = np.mean(all_y_vals)
    x_0 = np.mean(all_x_vals)

    distances = np.sqrt((np.array(all_x_vals) - x_0)**2 + (np.array(all_y_vals) - y_0)**2)
    max_distance = np.ceil(np.max(distances))

    # distances = distances/max_distance  # normalize distances to [0, 1]

    number_of_bins = 7

    bins = np.linspace(0, 1.05*np.max(distances), number_of_bins + 1)

    def make_stacked(fractions, distances, bins):
        stacked = [[] for _ in range(len(bins) - 1)]
        for dist, val in zip(distances, fractions):
            idx = np.searchsorted(bins, dist, side='right') - 1
            if 0 <= idx < len(stacked):
                stacked[idx].append(val)
        return stacked

    stacked_sf = make_stacked(all_surface_fractions, distances, bins)
    stacked_af = make_stacked(all_area_fractions, distances, bins)

    labels = [f'{int(bins[i])}-{int(bins[i+1])}' for i in range(len(bins)-1)]


    # Surface fraction stacked histogram
    xlims1 = (-0.01, 1.01)
    cmap = plt.cm.viridis_r
    if plot_mean_std:
        # compute the mean surface fraction for total distribution
        mean_surface_fraction = np.mean(all_surface_fractions)
        axs[iter, 0].axvline(x=mean_surface_fraction, color='red', linestyle='--', label=rf'$\mu_{{\gamma}}$ = {mean_surface_fraction:.2f}')
        standarddev_surface_fraction = np.std(all_surface_fractions)
        axs[iter, 0].axvline(x=mean_surface_fraction + standarddev_surface_fraction, color='blue', linestyle='--', label=rf'$\sigma_{{\gamma}}$ = {standarddev_surface_fraction:.2f}')
        axs[iter, 0].axvline(x=mean_surface_fraction - standarddev_surface_fraction, color='blue', linestyle='--')

     # Surface fraction stacked histogram

    axs[iter, 0].hist(stacked_sf, bins=np.arange(xlims1[0], xlims1[1], 0.01), stacked=True, label=labels, edgecolor='none', color=cmap(np.linspace(0, 1, len(labels))), density=False)
    if iter == 0:
        axs[iter, 0].set_title('PDF of Surface Fractions')
    if iter == len(framework_heads)-1:
        axs[iter, 0].set_xlabel('Surface Fraction (γ)')
    axs[iter, 0].set_ylabel('Number of Cells')
    axs[iter, 0].legend(title='Distance', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
    axs[iter, 0].set_xlim(xlims1)
    axs[iter, 0].set_yscale('log')

    # Surface fraction stacked histogram
    # xlims1 = (0.7, 1.001)
    cmap = plt.cm.viridis_r
    axs[iter, 1].hist(stacked_sf, bins=np.arange(xlims1[0], xlims1[1], 0.01), stacked=True, label=labels, edgecolor='none', color=cmap(np.linspace(0, 1, len(labels))), density=False, cumulative=True)
    if iter == len(framework_heads)-1:
        axs[iter, 1].set_xlabel('Surface Fraction (γ)')
    if iter == 0:
        axs[iter, 1].set_title('CDF of Surface Fractions', fontsize=10)
    axs[iter, 1].legend(title='Distance', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
    axs[iter, 1].set_xlim(xlims1)


    xlims2 = (-0.01, 1.01)
    cmap = plt.cm.inferno_r
    
    if plot_mean_std:
        # compute the mean area fraction for total distribution
        mean_area_fraction = np.mean(all_area_fractions)
        standardev_area_fraction = np.std(all_area_fractions)
        axs[iter, 2].axvline(x=mean_area_fraction, color='red', linestyle='--', label=rf'$\mu_{{\beta}}$ = {mean_area_fraction:.2f}')
        axs[iter, 2].axvline(x=mean_area_fraction + standardev_area_fraction, color='blue', linestyle='--', label=rf'$\sigma_{{\beta}}$ = {standardev_area_fraction:.2f}')
        axs[iter, 2].axvline(x=mean_area_fraction - standardev_area_fraction, color='blue', linestyle='--')

    # Area fraction stacked histogram
    axs[iter, 2].hist(stacked_af, bins=np.arange(xlims2[0], xlims2[1], 0.01), stacked=True, label=labels, edgecolor='none', color=cmap(np.linspace(0, 1, len(labels))), density=False)
    if iter == 0:
        axs[iter, 2].set_title('PDF of Area Fractions')
    if iter == len(framework_heads)-1:
        axs[iter, 2].set_xlabel('Area Fraction (β)')
    axs[iter, 2].legend(title='Distance', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
    axs[iter, 2].set_xlim(xlims2)

    # xlims2 = (0.0, 1.01)
    cmap = plt.cm.inferno_r
    # Area fraction stacked histogram
    axs[iter, 3].hist(stacked_af, bins=np.arange(xlims2[0], xlims2[1], 0.01), stacked=True, label=labels, edgecolor='none', color=cmap(np.linspace(0, 1, len(labels))), density=False, cumulative=True)
    
    if iter == 0:
        axs[iter, 3].set_title('CDF of Area Fractions', fontsize=10)
    if iter == len(framework_heads)-1:
        axs[iter, 3].set_xlabel('Area Fraction (β)')
    axs[iter, 3].legend(title='Distance', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
    axs[iter, 3].set_xlim(xlims2)

    # make a title for each row on the leftmost plot
    axs[iter, 0].set_ylabel(framework_head, fontsize=12, fontweight='bold', rotation=90, labelpad=40)

plt.tight_layout()
plt.show()

# save the figure
fname = f'MonolayerGrowth_Surface_Area_Fractions_Comparison'
fig.savefig('1000_Data/' + fname + '.png', dpi=800, bbox_inches='tight')




