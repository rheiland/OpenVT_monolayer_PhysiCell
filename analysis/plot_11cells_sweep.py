# Examples (run from directory containing the .mat files):
#   - plot (time,xpos) of the right-most cell in the 11-cell mechanics test
#
#  At t=0, an immovable "wall" cell is at x=0 and the remaining 10 cells overlap by 1 radius width
#  so that the right-most cell (ID=10) is at x=50. We let the cells relax (adhesion=0; repulsion varies)
#  and plot the curve with the right-most cell reaches x=90.
#

import sys
import glob
import os
import xml.etree.ElementTree as ET
import math
from pathlib import Path

from pyMCDS import pyMCDS
try:
  import matplotlib
  from matplotlib import gridspec
  import matplotlib.colors as mplc
  from matplotlib.patches import Circle, Ellipse, Rectangle
  from matplotlib.collections import PatchCollection
except:
  print("\n---Error: cannot import matplotlib")
  print("---Try: python -m pip install matplotlib")
#  print("---Consider installing Anaconda's Python 3 distribution.\n")
  raise
try:
  import numpy as np  # if mpl was installed, numpy should have been too.
except:
  print("\n---Error: cannot import numpy")
  print("---Try: python -m pip install numpy\n")
  raise
from collections import deque
try:
  # apparently we need mpl's Qt backend to do keypresses 
  matplotlib.use("Qt5Agg")
#   matplotlib.use("TkAgg")
  import matplotlib.pyplot as plt
except:
  print("\n---Error: cannot use matplotlib's TkAgg backend")
#  print("Consider installing Anaconda's Python 3 distribution.")
  raise

# current_idx = 0
nargs = len(sys.argv)-1
print("# args=",nargs)
max_idx = 1
if nargs > 0:
    max_idx = int(sys.argv[1])
print("max_idx= ",max_idx)

#for idx in range(len(sys.argv)):
use_defaults = True
show_nucleus = 0
current_idx = 0
axes_min = 0.0
axes_max = 1000  

current_idx = 0
print("current_idx=",current_idx)

#d={}   # dictionary to hold all (x,y) positions of cells

""" 
--- for example ---
In [141]: d['cell1599'][0:3]
Out[141]: 
array([[ 4900.  ,  4900.  ],
       [ 4934.17,  4487.91],
       [ 4960.75,  4148.02]])
"""

# fig = plt.figure(figsize=(7,5))
fig = plt.figure(figsize=(5,5))  # square
ax0 = fig.gca()


#-----------------------------------------------------
def get_cells_xpos(out_dir):
    global current_idx, axes_max,cax2,ax0,tvals,xpos
    # global current_idx, axes_max,cax2,ax0

    frame = current_idx 

    xml_file_root = "output%08d.xml" % frame
    # print("------ plot_cells_xpos():  current_idx= ",current_idx)
    # print("xml_file_root = ",xml_file_root)
    # xml_file = os.path.join('.', xml_file_root)
    # xml_file = os.path.join('output_11cells_symm', xml_file_root)
    # xml_file = os.path.join('output_11cells_symm_repuls_5', xml_file_root)
    # xml_file = os.path.join('output_11cells_symm_repuls_10', xml_file_root)
    # xml_file = os.path.join('output_optimize_11cells', 'out_run_001', xml_file_root)
    xml_file = os.path.join('output_optimize_11cells', out_dir, xml_file_root)
    # print("xml_file= ",xml_file)

    if not Path(xml_file).is_file():
        print("ERROR: file not found",xml_file)
        # return -1

    # mcds = pyMCDS(xml_file_root, microenv=False, graph=False, verbose=True)
    mcds = pyMCDS(xml_file, microenv=False, graph=False, verbose=False)
    total_min = mcds.get_time()  # warning: can return float that's epsilon from integer value
    # print("------ plot_cells_xpos():  total_min= ",total_min)
    try:
        df_all_cells = mcds.get_cell_df()
    except:
        print("plot_cells_xpos(): error performing mcds.get_cell_df()")
        sys.exit()
        # return
        
    xvals = df_all_cells['position_x']
    # print("type(xvals)= ",type(xvals))  # <class 'pandas.core.series.Series'>
    # print("xvals= ",xvals)

    # yvals = df_cells['position_y']

    xpos.append(xvals/10)   # divide to get units of cell diam
    # print("xpos= ",xpos)
            

    axes_min = mcds.get_mesh()[0][0][0][0]
    axes_max = mcds.get_mesh()[0][0][-1][0]

    # title_str = '11 horizontal cells mechanics test (PhysiCell)'
    # ax0.set_title(title_str, fontsize=12)

    # return 0

print("\nNOTE: click in plot window to give it focus before using keys.")


#----------------------------------------
out_dirs = ["out_run_000"]
out_dirs = ["out_run_001"]
out_dirs = ["out_run_000","out_run_001"]
out_dirs = ["out_run_000","out_run_001","out_run_002","out_run_003"]
out_dirs = []
for idx in range(9):
    out_dirs.append(f'out_run_{idx:03d}')

xpos = []
for out_dir in out_dirs:
    print("----------- ",out_dir)

    config_file = os.path.join('output_optimize_11cells', out_dir, "config.xml")
    tree = ET.parse(config_file)
    xml_root = tree.getroot()
    repulsion = xml_root.find('.//' + 'cell_cell_repulsion_strength').text
    adhesion = xml_root.find('.//' + 'cell_cell_adhesion_strength').text
    print("repulsion = ",repulsion)
    print("adhesion = ",adhesion)

    tvals = []
    # tvals.clear()
    xpos.clear()
    # max_idx = 577
    # max_idx = 5  # debugging
    # max_idx = 106
    # print("\n------------------- 1st loop: get tvals ")
    for idx in range(0,max_idx):
        xml_file_root = "output%08d.xml" % idx
        # print("---------- xml_file_root = ",xml_file_root)
        # xml_file = os.path.join('output_optimize_11cells', 'out_run_001', xml_file_root)
        xml_file = os.path.join('output_optimize_11cells', out_dir, xml_file_root)
        # print("---------- xml_file= ",xml_file)
        # mcds = pyMCDS(xml_file_root, microenv=False, graph=False, verbose=True)
        try:
            mcds = pyMCDS(xml_file, microenv=False, graph=False, verbose=False)
        except:
            break
        total_min = mcds.get_time()  # warning: can return float that's epsilon from integer value
        # print("total_min= ",total_min)
        tvals += [total_min]

    # print("\n------------------- 2nd loop: get xpos ")
    for idx in range(0,max_idx):
        current_idx = idx
        get_cells_xpos(out_dir)
        # ret_val = get_cells_xpos(out_dir)
        # if ret_val < 0:
        #    break

    # plt.plot(tvals,xpos,'o-', markersize=4)

    # Scale so a time unit=1 represents 90% relaxation. This will become the cell cycle duration for the monolayer.
    # t_90pct = 620.0
    #t_90pct = 443.0
    #t_90pct = 177.3
    t_90pct = 88.7
    t_90pct = 999

    # t_90_file = os.path.join('output_optimize_11cells', 'out_run_001',"time_90pct.txt")
    t_90_file = os.path.join('output_optimize_11cells', out_dir, "time_90pct.txt")
    try:
        with open(t_90_file, 'r') as file:
            t_90pct = float(file.readline())
        print("t_90pct = ", t_90pct)
    except:
        print(f'--- Error: unable to read {t_90_file}')
        pass

    tvals_ = np.array(tvals)   
    xpos_ = np.array(xpos)   
    len_tvals = len(tvals)
    # print("len(tvals)=", len_tvals)
    # print("len(xpos)=",len(xpos))
    # tv = tvals_/t_90pct
    # xv_start = xpos_[:,0]
    # print("xv_start =",xv_start)
    xv_end = xpos_[:,10]
    # print("xv_end =",xv_end)
    tissue_width = xpos_[:,10] - xpos_[:,0]
    if False:   # rwh
        with open("pc_plot_11cells.csv", 'w') as f:
            for idx in range(len(tv)):
                f.write(f'{tv[idx]},{tissue_width[idx]}\n')
        f.close()

    # plt.plot(tvals_/t_90pct, tissue_width,'-', markersize=4)   # only plot the "last" curve (right-most cell)
    # if out_dir == "out_run_000":
    #     plt.plot(tvals_/t_90pct, tissue_width,'--',color='g')   # only plot the "last" curve (right-most cell)
    # else:
    #     plt.plot(tvals_/t_90pct, tissue_width,'r')   # only plot the "last" curve (right-most cell)

    # print(np.array2string(tvals_/t_90pct, separator=','))
    # print(np.array2string(tissue_width, separator=','))

    # plt.plot(tvals_/t_90pct, tissue_width,'.',markersize=4, label=out_dir)
    # plt.plot(tvals_/t_90pct, tissue_width, label=out_dir)
    t_calib = tvals_/t_90pct
    # plt.plot(tvals_/t_90pct, tissue_width, label=f'{adhesion},{repulsion}')
    plt.plot(t_calib, tissue_width, label=f'{adhesion},{repulsion}')
    xtext = (tvals_/t_90pct)[:-1][:-1]
    ytext = tissue_width[:-1][:-1]
    xtext = t_calib[-1]
    ytext = tissue_width[-1]
    print("t_calib=",t_calib)
    print("len(t_calib)=",len(t_calib))
    print("tissue_width=",tissue_width)
    plt.text(xtext, ytext, f'{adhesion},{repulsion}')
    # plt.annotate(f'{adhesion},{repulsion}', xy=(3, 6), xytext=(2.5, 5.5)
    # plt.show()

# plt.plot(tvals_/t_90pct, np.sin(tvals_) )
#------------------------
# ax0.set_xlim(0, 5)
ax0.set_xlim(0, 10)
ax0.set_ylim(5, 10)

# draw horiz and vertical dashed lines for rightmost cell reaching 90% relaxation width
plt.plot([0,10],[9,9],'--k')
plt.plot([1,1],[0,10],'--k')  # if scaled to "CD"

x, y = np.loadtxt('analysis/relaxation_exact.csv', skiprows=1, delimiter=',',  unpack=True)
plt.plot(x, y,color='blue',linestyle='--',label='analytic')

ax0.set_xlabel("Time (relative to 90% width)", fontsize=14)
ax0.set_ylabel("Tissue width (CD)", fontsize=14)
title_str = 'PhysiCell: relaxation of 11 cells'
ax0.set_title(title_str, fontsize=12)

plt.legend()

plt.show()
