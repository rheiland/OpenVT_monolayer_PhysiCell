# cf. montage_monolayer.py

import os
#import subprocess

output_dirs = []

montage_cmd = "montage -geometry +0+0 -tile 5x1 "

# NEW: 4/23/26
gamma_vals= [0.0, 0.696, 0.819, 0.895, 0.903]
beta_vals= [ 0.995]

print("type(gamma_vals)=", type(gamma_vals))
print(gamma_vals)
gamma_r = list(reversed(gamma_vals))
print("type(gamma_r)=", type(gamma_r))
print(gamma_r)
for beta in beta_vals:
    # for gamma in gamma_r:
    for gamma in gamma_vals:
        output_dir = "out_cell_area_b" + str(beta) + "_g" + str(gamma)
        # cmd = "python ../beta/plot_cell_scalars-2.py -o " + output_dir + " -s beta_or_gamma -f -1 -a -x0 -700 -x1 700 -y0 -700 -y1 700"
        cmd = "python ../beta/plot_cell_scalars_4states.py -o " + output_dir + " -s beta_or_gamma -f -1 -a -x0 -750 -x1 750 -y0 -750 -y1 750"
        os.system(cmd)

        montage_cmd += output_dir + "/keep.png "

# montage -geometry +0+0 -tile 5x5 run01/keep.png run02/keep.png run03/keep.png run04/keep.png run05/keep.png row1.png
montage_cmd += " montage_1x5.png"
print("montage_cmd= ",montage_cmd)