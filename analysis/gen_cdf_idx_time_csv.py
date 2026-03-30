
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
max_dirs = int(sys.argv[1])
tau = float(sys.argv[2])   # = 5*T (T=time to reach 90% width of 11-cells relaxation model); e.g., 443.5
print("doing max_dirs= ",max_dirs)
print("tau= ",tau)

# cell_radius = 5.0

# t=[]
# tumor_diam=[]
# fig, ax = plt.subplots()

# for irun in range(100):
for irun in range(max_dirs):
# for irun in [22,31,39,40,50,63,83]:
# for irun in [14]:
    # out_dir = f'out_cells1000_{irun}'
    out_dir = f'bg00_cells1000_{irun}'
    print("out_dir= ",out_dir)
    # print("xml_file= ",xml_file)

    run_dir = "run" + str(irun+1)
    if (not os.path.exists(run_dir)):
        print("--- mkdir ", run_dir)
        os.makedirs(run_dir)

    xml_pattern = out_dir + "/" + "output*.xml"
    xml_files = glob.glob(xml_pattern)
    xml_files.sort()
    # last_file = xml_files[-1]
    # print("last_file= ",last_file)

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

        # print("cells_x/cell_radius= ",cells_x/cell_radius)

        # current_time = mcds.get_time()
        t_normd.append(mcds.get_time() / tau)
        # print('time (min)= ', current_time )
        # print('time (hr)= ', current_time/60. )
        # print('time (day)= ', current_time/1440. )

        # print("# cells= ",cells_x_calibrated.shape[0])
        # diam = cells_x.max() - cells_x.min()
        # print("monolayer diam= ",diam)

        # t.append(current_time/1440.)  # to get days
        # t.append(current_time/88.7)  # recall the 88.7 mins (the 90% width of 11 cells) = 1 T unit 
        # tumor_diam.append(diam/5.0)      # calibrate space units by dividing by cell radius
        # sub_intern.append(sintern)
        # sub_conc.append(sconc)

    # file_out = f'{out_dir}/cell_data_no_inhibition_{irun}.csv'
    # file_out = f'cell_data_no_inhibition_{irun}.csv'
    # file_out = f'timeline.txt'
    file_out  = os.path.join(run_dir,'timeline.txt')
    print("--> ",file_out)
    with open(file_out, "w", newline="") as file:
    # with open(os.path.join(run_dir,file_out), "w", newline="") as file:
        writer = csv.writer(file, delimiter=' ')

        writer.writerow(['index','t'])
        # writer.writerow(['x_pos','y_pos','radius_i','f_i','a_i'])
        for jdx in range(len(t_normd)):
            formatted_row = (jdx, f'{t_normd[jdx]:.6f}')
            # writer.writerow([jdx,t_normd[jdx]])
            writer.writerow(formatted_row)
