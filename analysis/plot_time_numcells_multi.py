# Copy this script and pyMCDS.py into the output directory (/output_unit_test1)
#
# e.g.:
# python analyze_monolayer.py 330 331
#
__author__ = "Randy Heiland"

import sys,pathlib
#import xml.etree.ElementTree as ET
#import math
# import scipy.io
# from pyMCDS import pyMCDS
from pyMCDS_cells import pyMCDS_cells
#import matplotlib
#import numpy as np
import matplotlib.pyplot as plt

# print(len(sys.argv))
# if (len(sys.argv) < 4):
#     print("--- No args provided, will try to use 'output' dir and 0th output file")  
#     out_dir = "output"
#     idx_min = 0
#     idx_max = 1
# else:
#     kdx = 1
#     out_dir = sys.argv[kdx]
#     print("out_dir=",out_dir)
#     kdx += 1
#     idx_min = int(sys.argv[kdx])
#     kdx += 1
#     idx_max = int(sys.argv[kdx])

# print('frame, field = ',frame_idx, field_index)

run_dirs = ["run1", "run2", "run3", "run4"] 
run_dirs = ["run1", "run2", "run3", "run4", "run5", "run6"] 
run_colors = ["red", "green", "black", "orange", "cyan", "purple"] 

#out_dir = "../PhysiCell/output_monolayer_pressure_set_behavior"
tumor_diam=[]

fig, ax = plt.subplots()

# ------- 1st plot all computed values (at every 10 hours)
hr_delta = 1
#hr_delta = 10
max_frame = 1000
# for out_dir in run_dirs for run_color in run_colors:
t=[]
tumor_numcells=[]
for out_dir,run_color in zip(run_dirs,run_colors):
  t.clear()
  tumor_numcells.clear()
  for idx in range(0,max_frame+1, hr_delta):
    xml_file = "output%08d.xml" % idx
    print("xml_file= ",xml_file)

    # mcds = pyMCDS(xml_file, '../PhysiCell/output_usecase0')   # reads BOTH cells and substrates info
    try:
        # mcds = pyMCDS(xml_file, out_dir)   # reads BOTH cells and substrates info
        mcds = pyMCDS_cells(xml_file, out_dir)   # reads BOTH cells and substrates info
    except:
        break
    # cell_ids = mcds.data['discrete_cells']['ID']
    cells_x = mcds.data['discrete_cells']['position_x']

    current_time = mcds.get_time()
    # print('time (min)= ', current_time )
    # print('time (hr)= ', current_time/60. )
    # print('time (day)= ', current_time/1440. )
    num_cells = cells_x.shape[0]
    print("# cells= ",cells_x.shape[0])
    # if num_cells > 10000:   # stop 1 output step short
    #    break
    # diam = cells_x.max() - cells_x.min()
    # print("monolayer diam= ",diam)

    # t.append(current_time/1440.)  # to get days
    t.append(current_time/88.7)  # recall the 88.7 mins (the 90% width of 11 cells) = 1 T unit 
    tumor_numcells.append(num_cells)      # calibrate space units by dividing by cell radius
    # sub_intern.append(sintern)
    # sub_conc.append(sconc)

  print("t=",t)
  print("tumor_numcells=",tumor_numcells)
  ax.plot(t, tumor_numcells,'k-', label=out_dir, color=run_color)
  ax.plot(t, tumor_numcells,'k.')
  print('--------------\n')


# ax.set_xlim(0, 110)
# ax.set_ylim(0, 190)

ax.legend(loc='lower right')

ax.set(xlabel='time (T)', ylabel='# cells',title="monolayer growth (PhysiCell)")
# ax.grid()
# fig.savefig("test.png")
plt.show()
