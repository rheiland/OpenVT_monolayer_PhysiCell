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

out_dir = "../PhysiCell/output_monolayer"
out_dir = "../PhysiCell/output_monolayer_final"
out_dir = "./output_monolayer2"
out_dir = "./output_monolayer_live_cycle"
out_dir = "./output_beta_default_30day"
out_dir = "./output_beta_14days"
out_dir = "./output_beta_14days_fixcycle_18hr"

out_dir = "./output_monolayer_live_cycle_stochastic"
out_dir = "./output_monolayer_cytomsep"
out_dir = "./output_monolayer"

#out_dir = "../PhysiCell/output_monolayer_pressure_set_behavior"
t=[]
tumor_diam=[]

fig, ax = plt.subplots()

# ------- 1st plot all computed values (at every 10 hours)
hr_delta = 20
hr_delta = 1
#hr_delta = 10
# for idx in range(0,648, hr_delta):
max_frame = 615
max_frame = 162
max_frame = 120
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
    print('time (hr)= ', current_time/60. )
    # print('time (day)= ', current_time/1440. )
    print("# cells= ",cells_x.shape[0])
    diam = cells_x.max() - cells_x.min()
    print("monolayer diam= ",diam)

    t.append(current_time/1440.)
    tumor_diam.append(diam)
    # sub_intern.append(sintern)
    # sub_conc.append(sconc)

ax.plot(t, tumor_diam,'k-')

t2 = []
tumor_diam2 = []
add_points_flag = True
add_points_flag = False
# ------- 2nd plot just the points from the table
drasdo_time=[336/24.,386/24.,408/24.,481/24.,506/24.,646/24.]
drasdo_diam=[1140,1400,1590,2040,2250,3040]
if add_points_flag:
#   for idx in [336, 386, 408, 481, 506, 646]:   # times (hours) from the table
#   for idx in [336, 386, 408, 481, 506, 614]:   # times (hours) from the table
  for idx in [84, 96, 102, 120, 126, 153]:   # approx if save interval =4 hrs
    xml_file = "output%08d.xml" % idx

    # mcds = pyMCDS(xml_file, '../PhysiCell/output_usecase0')   # reads BOTH cells and substrates info
    # mcds = pyMCDS(xml_file, out_dir)   # reads BOTH cells and substrates info
    mcds = pyMCDS_cells(xml_file, out_dir)   # reads BOTH cells and substrates info
    # cell_ids = mcds.data['discrete_cells']['ID']
    cells_x = mcds.data['discrete_cells']['position_x']

    current_time = mcds.get_time()
    # print('time (min)= ', current_time )
    print('time (hr)= ', current_time/60. )
    # print('time (day)= ', current_time/1440. )
    print("# cells= ",cells_x.shape[0])
    diam = cells_x.max() - cells_x.min()
    print("monolayer diam= ",diam)

    t2.append(current_time/1440.)
    tumor_diam2.append(diam)

if add_points_flag:
  ax.plot(t2, tumor_diam2,'ko')

ax.plot(drasdo_time, drasdo_diam,'ro')

ax.set(xlabel='t (day)', ylabel='diameter (micron)',title="monolayer growth")
# ax.grid()
# fig.savefig("test.png")
plt.show()
