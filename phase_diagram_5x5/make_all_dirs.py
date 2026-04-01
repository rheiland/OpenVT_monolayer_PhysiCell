import os
import sys
import subprocess

# NEW: 3/27/26
# (base) M1P~/git/PhysiCell_monolayer/PhysiCell_mech_grid_xml$ python beta/chosen_CDF_pct.py 100 cdf_1000cells_linear_growth  bg00_cells1000_  f_i
# -->
# x values for %s: 0.00000, 0.72514, 0.83539, 0.90324, 0.91286
gamma_vals = [0.0, 0.725, 0.835, 0.903, 0.913]
#gamma_vals = [ 0.903]

# and:
# (base) M1P~/git/PhysiCell_monolayer/PhysiCell_mech_grid_xml$ python beta/chosen_CDF_pct.py 100 cdf_1000cells_linear_growth  bg00_cells1000_  a_i
# -->
# x values for %s: 0.51875, 0.98036, 0.98991, 0.99496, 0.99631
beta_vals  = [0.0, 0.980, 0.990, 0.995, 0.996]   # use 0.0 instead of 0.xxx
#beta_vals  = [ 0.996]   # use 0.0 instead of 0.xxx

output_dirs = []
for beta in beta_vals:
    for gamma in gamma_vals:
        folder_name = "out_cell_area_b" + str(beta) + "_g" + str(gamma)
        output_dirs.append(folder_name)
        if (not os.path.exists(folder_name)):
            print("--- mkdir ", folder_name)
            os.makedirs(folder_name)

            cmd = f"cp out_cell_area_b0.996_g0.913/* {folder_name}/."
            print("----- cmd = ",cmd)
            os.system(cmd)   # <------ Execute the simulation

