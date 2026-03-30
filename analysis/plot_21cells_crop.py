# Examples (run from directory containing the .mat files):
#   - plot (time,xpos) of the right-most cell in the 21-cell mechanics test
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
# print("# args=",len(sys.argv)-1)
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

tvals = []
xpos = []
#-----------------------------------------------------
def get_cells_xpos():
    global current_idx, axes_max,cax2,ax0,tvals,xpos

    frame = current_idx 

    xml_file_root = "output%08d.xml" % frame
    # print("plot_cell_scalar():  current_idx= ",current_idx)
    # print("xml_file_root = ",xml_file_root)
    # xml_file = os.path.join('.', xml_file_root)

    # xml_file = os.path.join('output_21cells_symm', xml_file_root)
    # xml_file = os.path.join('output_21cells_symm2', xml_file_root)
    xml_file = os.path.join('output_21cells', xml_file_root)

    if not Path(xml_file).is_file():
        print("ERROR: file not found",xml_file)
        return

    mcds = pyMCDS(xml_file, microenv=False, graph=False, verbose=False)
    total_min = mcds.get_time()  # warning: can return float that's epsilon from integer value
    try:
        df_all_cells = mcds.get_cell_df()
    except:
        print("vis_tab.py: plot_cell_scalar(): error performing mcds.get_cell_df()")
        return
        
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

print("\nNOTE: click in plot window to give it focus before using keys.")

# max_idx = 5  # debugging
# max_idx = 577
# max_idx = 1440
for idx in range(0,max_idx):
    xml_file_root = "output%08d.xml" % idx
    # print("xml_file_root = ",xml_file_root)
    # xml_file = os.path.join('.', xml_file_root)

    # xml_file = os.path.join('output_21cells_symm', xml_file_root)
    # xml_file = os.path.join('output_21cells_symm2', xml_file_root)
    xml_file = os.path.join('output_21cells', xml_file_root)
    mcds = pyMCDS(xml_file, microenv=False, graph=False, verbose=False)
    total_min = mcds.get_time()  # warning: can return float that's epsilon from integer value
    # print("total_min= ",total_min)
    tvals += [total_min]

for idx in range(0,max_idx):
    current_idx = idx
    get_cells_xpos()

# plt.plot(tvals,xpos,'o-', markersize=4)

# Scale so a time unit=1 represents 90% relaxation. This will become the cell cycle duration for the monolayer.
# t_90pct = 620.0
# t_90pct = 443.0
t_90pct = 88.7
tvals_ = np.array(tvals)   
xpos_ = np.array(xpos)   
len_tvals = len(tvals)
# print("len(tvals)=", len_tvals)
# print("len(xpos)=",len(xpos))
# print("tvals_ =",tvals_)
# print("xpos_ =",xpos_)
# plt.plot(tvals_/t_90pct, xpos,'-', markersize=4)
# plt.plot(tvals_/t_90pct, xpos,'-o', markersize=4)
tv = tvals_/t_90pct
# xv = xpos_[:,10]
# xv = xpos_[:,20]
tissue_width = xpos_[:,20] - xpos_[:,0]
csv_file = "pc_plot_21cells_width.csv"
with open(csv_file, 'w') as f:
    # for idx in range(len(tv)):
    for idx in range(len(tissue_width)):
        # f.write(f'{tv[idx]},{xv[idx]}\n')
        f.write(f'{tv[idx]},{tissue_width[idx]}\n')
f.close()
print("----> ",csv_file)
# plt.plot(tvals_/t_90pct, xpos_[:,20],'-', markersize=4)   # only plot the "last" curve (right-most cell)
# plt.plot(tv,xv,'-', markersize=4)   # only plot the "last" curve (right-most cell)
plt.plot(tv,tissue_width,'-', markersize=4)   # only plot the "last" curve (right-most cell)

# ax0.set_xlim(0, 16)
ax0.set_xlim(0, 10)
ax0.set_ylim(14, 20)

# draw horiz and vertical dashed lines for rightmost cell reaching 90% relaxation width
# plt.plot([0,5],[9,9],'--k')
# plt.plot([1,1],[0,10],'--k')  # if scaled to "CD"

# ax0.set_xlabel("Time (min)", fontsize=14)
ax0.set_xlabel("Time (relative to 90% width)", fontsize=14)

# ax0.set_ylabel("Cell center (microns)", fontsize=14)
# ax0.set_ylabel("Position (CD)", fontsize=14)
ax0.set_ylabel("Tissue width (CD)", fontsize=14)

# title_str = '11 horizontal cells mechanics test (PhysiCell)'
# title_str = '11 horiz cells mechanics test (PhysiCell)'
title_str = 'PhysiCell: relaxation test (11+10 cells)'
ax0.set_title(title_str, fontsize=12)

# keep last plot displayed
#plt.ioff()
# ax0.set_aspect('equal')
plt.show()
