
# generate .csv data CDF plots (100 runs of 1000 cells)
#  "index,t"  where time, t, is normalized by 5T 
#             (T= time to reach 90% width of 11-cells relaxation))

import sys
import os
import pathlib
import csv
import glob
from pyMCDS_cells import pyMCDS_cells
# import matplotlib.pyplot as plt

print("sys.argv=",sys.argv)
out_dir = sys.argv[1]
tau = float(sys.argv[2])   # = 5*T (T=time to reach 90% width of 11-cells relaxation model); e.g., 443.5
print("tau= ",tau)

# cell_radius = 5.0

# t=[]
# tumor_diam=[]
# fig, ax = plt.subplots()

for irun in range(1):
    # run_dir = "run" + str(irun+1)
    run_dir = out_dir + "_time_series"
    if (not os.path.exists(run_dir)):
        print("--- mkdir ", run_dir)
        os.makedirs(run_dir)

    xml_pattern = out_dir + "/" + "output*.xml"
    xml_files = glob.glob(xml_pattern)
    xml_files.sort()
    last_xml_file = xml_files[-1]
    print("last_xml_file= ",last_xml_file)

    t_normd = []
    idx = 0
    for xml_file in xml_files:
        try:
            # mcds = pyMCDS(xml_file, out_dir)   # reads BOTH cells and substrates info
            mcds = pyMCDS_cells(os.path.basename(xml_file), out_dir)   # reads BOTH cells and substrates info
            # mcds = pyMCDS(xml_file_root, self.output_dir, microenv=False, graph=False, verbose=False)
            # df_cells = get_mcds_cells_df(mcds)
            # df_all_cells = mcds.get_cell_df()
        except:
            print(f"pyMCDS_cells error reading {out_dir}/{xml_file}")
            exit
        # cell_ids = mcds.data['discrete_cells']['ID']
        # print(mcds.data['discrete_cells'].keys())
        # x_pos = mcds.data['discrete_cells']['position_x']
        # y_pos = mcds.data['discrete_cells']['position_y']
        # radius_i = mcds.data['discrete_cells']['radius']
        # f_i = mcds.data['discrete_cells']['f_i']
        # a_i = mcds.data['discrete_cells']['a_i']

        t_normd.append(mcds.get_time() / tau)

    # file_out  = os.path.join(run_dir,f'data{idx:06}.csv')
    file_out  = os.path.join(run_dir,'timeline.txt')
    print("--> ",file_out)
    with open(file_out, "w", newline="") as file:
    # with open(os.path.join(run_dir,file_out), "w", newline="") as file:
        writer = csv.writer(file, delimiter=',')

        writer.writerow(['index','t'])
        for jdx in range(len(t_normd)):
            formatted_row = (jdx, f'{t_normd[jdx]:.6f}')
            # writer.writerow([jdx,t_normd[jdx]])
            writer.writerow(formatted_row)
