# cf. montage_monolayer.py

import os
#import subprocess

# output_dirs = ['out_b0.5_g0.5', 'out_b0.5_g0.6', 'out_b0.5_g0.7', 'out_b0.5_g0.8', 'out_b0.5_g0.9', 'out_b0.6_g0.5', 'out_b0.6_g0.6', 'out_b0.6_g0.7', 'out_b0.6_g0.8', 'out_b0.6_g0.9', 'out_b0.7_g0.5', 'out_b0.7_g0.6', 'out_b0.7_g0.7', 'out_b0.7_g0.8', 'out_b0.7_g0.9', 'out_b0.8_g0.5', 'out_b0.8_g0.6', 'out_b0.8_g0.7', 'out_b0.8_g0.8', 'out_b0.8_g0.9', 'out_b0.9_g0.5', 'out_b0.9_g0.6', 'out_b0.9_g0.7', 'out_b0.9_g0.8', 'out_b0.9_g0.9']
# output_dirs = ['out_b0.5_g0.5', 'out_b0.5_g0.6']
output_dirs = []

# montage_cmd = "montage -geometry +0+0 -tile 9x9 "
montage_cmd = "montage -geometry +0+0 -tile 5x5 "

# Roman had beta in X (columns), gamma in Y (rows)
# gamma_vals = [0.5, 0.6, 0.7, 0.8, 0.9]
# gamma_vals = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
# gamma_vals = [0.00000, 0.68266, 0.77982, 0.81597, 0.86907]
# gamma_vals = [0.00000, 0.68266, 0.77982, 0.81597, 0.9081]
# gamma_vals = [0.00000, 0.68266, 0.77982, 0.9081, 0.93674]
# gamma_vals = [0.0, 0.68266, 0.77982, 0.9081, 0.93674]

# beta_vals  = [0.0, 0.978, 0.988, 0.994, 0.996]   # use 0.0 instead of 0.499 ??
# gamma_vals = [0.0, 0.711, 0.826, 0.903, 0.911]
# gamma_vals = [0.0, 0.711, 0.826]

# NEW: 4/7/26
gamma_vals= [0.0, 0.786, 0.891, 0.931, 0.948]
beta_vals= [0.0, 0.988, 0.995, 0.998, 0.999]

    # for gamma in [0.9081]:
# gamma_vals = [ 0.86907]
print("type(gamma_vals)=", type(gamma_vals))
print(gamma_vals)
gamma_r = list(reversed(gamma_vals))
print("type(gamma_r)=", type(gamma_r))
print(gamma_r)
for gamma in gamma_r:
    # for beta in [0.5, 0.6, 0.7, 0.8, 0.9]:
    # for beta in [0.09078, 0.76600, 0.94701, 0.97044, 0.98660 ]:
    # for beta in [0.09078, 0.766, 0.94701, 0.97044, 0.9866 ]:
    for beta in beta_vals:
    # for beta in [ 0.98660 ]:
        output_dir = "out_cell_area_b" + str(beta) + "_g" + str(gamma)
        # cmd = "python ../beta/plot_cell_scalars-2.py -o " + output_dir + " -s beta_or_gamma -f -1 -a -x0 -700 -x1 700 -y0 -700 -y1 700"
        cmd = "python ../beta/plot_cell_scalars_4states.py -o " + output_dir + " -s beta_or_gamma -f -1 -a -x0 -750 -x1 750 -y0 -750 -y1 750"
        os.system(cmd)

        montage_cmd += output_dir + "/keep.png "

# montage -geometry +0+0 -tile 5x5 run01/keep.png run02/keep.png run03/keep.png run04/keep.png run05/keep.png row1.png
montage_cmd += " montage.png"
print("montage_cmd= ",montage_cmd)