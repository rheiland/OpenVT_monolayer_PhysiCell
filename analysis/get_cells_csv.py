
import sys
import os
import pathlib
import csv
import glob
from pyMCDS_cells import pyMCDS_cells
import matplotlib.pyplot as plt

# out_dir = "./output_monolayer"
# out_dir = "./run02_nbrs"

print("sys.argv=",sys.argv)
out_dir = sys.argv[1]

# Jan_m1.log:custom_function: 120: cell ID= 10, volume= 523.6, radius= 5

# total volume
#(base) M1P~/git/monolayer$ grep "<total" run01/config.xml
#                    <total units="micron^3">523.6</total>
#(base) M1P~/git/monolayer$ grep "<total" run25/config.xml
#                    <total units="micron^3">523.6</total>
cell_radius = 5.0

# t=[]
# tumor_diam=[]
# fig, ax = plt.subplots()
# ------- 1st plot all computed values (at every 10 hours)
hr_delta = 1
# for idx in range(1,2, hr_delta):
idx = 5
# xml_file = "output%08d.xml" % idx
print("out_dir= ",out_dir)
# print("xml_file= ",xml_file)

xml_pattern = out_dir + "/" + "output*.xml"
xml_files = glob.glob(xml_pattern)
xml_files.sort()
last_file = xml_files[-1]
print("last_file= ",last_file)

try:
    # mcds = pyMCDS(xml_file, out_dir)   # reads BOTH cells and substrates info
    mcds = pyMCDS_cells(os.path.basename(last_file), out_dir)   # reads BOTH cells and substrates info
    # mcds = pyMCDS(xml_file_root, self.output_dir, microenv=False, graph=False, verbose=False)
    # df_cells = get_mcds_cells_df(mcds)
    # df_all_cells = mcds.get_cell_df()
except:
    print("pyMCDS_cells error reading ",out_dir,last_file)
    exit
# cell_ids = mcds.data['discrete_cells']['ID']
# print(mcds.data['discrete_cells'].keys())
cells_x_calibrated = mcds.data['discrete_cells']['position_x'] / cell_radius
cells_y_calibrated = mcds.data['discrete_cells']['position_y'] / cell_radius
beta_or_gamma = mcds.data['discrete_cells']['beta_or_gamma']   # 0: growing; >0: arrested
num_nbrs = mcds.data['discrete_cells']['num_nbrs']
# print("cells_x_calibrated= ",cells_x_calibrated)
# print("cells_y_calibrated= ",cells_y_calibrated)

# print("cells_x/cell_radius= ",cells_x/cell_radius)

current_time = mcds.get_time()
# print('time (min)= ', current_time )
# print('time (hr)= ', current_time/60. )
# print('time (day)= ', current_time/1440. )
print("# cells= ",cells_x_calibrated.shape[0])
# diam = cells_x.max() - cells_x.min()
# print("monolayer diam= ",diam)

# t.append(current_time/1440.)  # to get days
# t.append(current_time/88.7)  # recall the 88.7 mins (the 90% width of 11 cells) = 1 T unit 
# tumor_diam.append(diam/5.0)      # calibrate space units by dividing by cell radius
# sub_intern.append(sintern)
# sub_conc.append(sconc)
file_out = "cells.csv"
print("--> ",file_out)
with open(file_out, "w", newline="") as file:
    writer = csv.writer(file)

# Write each row: x,y,g,n  (where g=growing (0/1), n=# of nbrs)
    for jdx in range(len(cells_x_calibrated)):
        if beta_or_gamma[jdx] == 0.0:   # growing, but confusingly, Roman maps 0=inhibited and 1=growing
            growing = 1
        else:
            growing = 0
        # num_nbrs = 6  # dummy for now
        # num_nbrs = 3  # dummy for now
        writer.writerow([cells_x_calibrated[jdx],cells_y_calibrated[jdx],growing,int(num_nbrs[jdx])])

