
# generate .csv data CDF plots (100 runs of 1000 cells) for Dom to comparing framework results

import sys
import os
import pathlib
import csv
import glob
from pyMCDS_cells import pyMCDS_cells
import matplotlib.pyplot as plt

# print("sys.argv=",sys.argv)
# out_dir = sys.argv[1]

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

out_dir = 'output_linear_growth'
if True:
# for irun in range(999):
# for irun in [22,31,39,40,50,63,83]:
# for irun in [14]:
    # out_dir = f'out_cells1000_{irun}'
    # out_dir = f'bg00_cells1000_{irun}'
    # print("out_dir= ",out_dir)
    # print("xml_file= ",xml_file)

    xml_pattern = out_dir + "/" + "output*.xml"
    xml_files = glob.glob(xml_pattern)
    xml_files.sort()
    last_file = xml_files[-1]
    # print("last_file= ",last_file)
    print("xml_files= ",xml_files)
    idx_time = 0
    for xml_file in xml_files:

        file_out = f'cell_data_no_inhibition_{idx_time}.csv'
        idx_time += 1
        print("--> ",file_out)
        with open(file_out, "w", newline="") as file:
            writer = csv.writer(file)

            try:
                # mcds = pyMCDS(xml_file, out_dir)   # reads BOTH cells and substrates info
                mcds = pyMCDS_cells(os.path.basename(xml_file), out_dir)   # reads BOTH cells and substrates info
                # mcds = pyMCDS(xml_file_root, self.output_dir, microenv=False, graph=False, verbose=False)
                # df_cells = get_mcds_cells_df(mcds)
                # df_all_cells = mcds.get_cell_df()
            except:
                print("pyMCDS_cells error reading ",out_dir,last_file)
                exit

            # current_time = mcds.get_time()

            # cell_ids = mcds.data['discrete_cells']['ID']
            # print(mcds.data['discrete_cells'].keys())
            x_pos = mcds.data['discrete_cells']['position_x']
            y_pos = mcds.data['discrete_cells']['position_y']
            radius_i = mcds.data['discrete_cells']['radius']
            f_i = mcds.data['discrete_cells']['f_i']
            a_i = mcds.data['discrete_cells']['a_i']

            # file_out = f'{out_dir}/cell_data_no_inhibition_{irun}.csv'

            # Write each row: x,y,g,n  (where g=growing (0/1), n=# of nbrs)
            writer.writerow(['x_pos','y_pos','radius_i','f_i','a_i'])
            for jdx in range(len(x_pos)):
                writer.writerow([x_pos[jdx],y_pos[jdx],radius_i[jdx],f_i[jdx],a_i[jdx]])
