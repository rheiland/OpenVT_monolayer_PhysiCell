
import sys
import os
import glob
from pyMCDS import pyMCDS
import numpy as np
import matplotlib.pyplot as plt

def fix_hist_step_vertical_line_at_end(ax):
    """
    Helper function to remove the vertical line artifact from a step-type histogram.
    """
    for poly in ax.findobj(plt.Polygon):
        # Retrieve the vertices and remove the last one, which is the final vertical line
        new_vertices = poly.get_xy()[:-1]
        poly.set_xy(new_vertices)

argc=len(sys.argv)
print("argc= ",argc)
if argc < 5:
    print("Missing <num_runs> <cdf_dir> <dir_prefix> <f_i or a_i>")
    print("e.g., 100 cdf_1000cells out_cells1000_  a_i")
    sys.exit()

print('argv[0]=',sys.argv[0])
# cell_scalar_name = sys.argv[1]
print('argv[1]=',sys.argv[1])
max_runs = int(sys.argv[1])
output_dir_base = sys.argv[2]
dir_prefix = sys.argv[3]
cell_scalar_name = sys.argv[4]

# max_runs = 2
# max_runs = 100
# max_runs = 10
# output_dir_base = "cdf_100cells"
# total_count = max_runs * 100
# suffix = "1000"
# suffix = "200"
# output_dir_base = "cdf_1000cells"
# output_dir_base = "cdf_200cells"
# output_dir_base = f"cdf_{suffix}cells"
# output_dir_base += "_slowgrowth_vol2494"
# total_count = max_runs * 1000   # turns out we don't need/want this!
# cell_scalar_name = "f_i"    # surface fraction overlap with nbrs  (~gamma)
# cell_scalar_name = "a_i"    # area fraction overlap with nbrs  (~beta)

# fig, ax = plt.subplots()
# fig = plt.figure(figsize=(7, 5))
fig = plt.figure(figsize=(5, 5))

for idx in range(max_runs):  # for each of the dirs/replicate runs
    xml_file_root = "output%08d.xml" % idx
    # output_dir = output_dir_base + f'/bg00_cells100_{idx}'
    # output_dir = output_dir_base + f'/bg00_cells100_{idx}'
    # if output_dir_base == "cdf_1000cells":
    #     output_dir = output_dir_base + f'/bg00_cells1000_{idx}'
    # elif output_dir_base == "cdf_200cells":
    #     output_dir = output_dir_base + f'/bg00_cells200_{idx}'
    # output_dir = output_dir_base + f'/bg00_cells{suffix}_{idx}'
    output_dir = f'{output_dir_base}/{dir_prefix}{idx}'
    print(idx,")  ------- output_dir= ",output_dir)
    xml_pattern = output_dir + "/" + "output*.xml"
    xml_files = glob.glob(xml_pattern)
    xml_files.sort()
    print("xml_files= ",xml_files)
    last_file = xml_files[-1]
    print("-------- last_file =",last_file )
    xml_file_root = os.path.basename(last_file)
    # print("xml_file_root =",xml_file_root )
    print
    mcds = pyMCDS(xml_file_root, output_dir, microenv=False, graph=False, verbose=False)
    # df_cells = self.get_mcds_cells_df(mcds)
    df_cells = mcds.get_cell_df()
    cell_scalar = df_cells[cell_scalar_name]
    # if cell_scalar.max() < 0.4:
    #     print("max = ",cell_scalar.max())
    if idx == 0:
        all_vals = cell_scalar
    else:
        # all_vals = np.concatenate((all_vals, cell_scalar), axis=0)
        all_vals = np.concatenate((all_vals, cell_scalar), axis=0)
    # print(cell_scalar)

# print("-------- last_file =",last_file )
# print("len(all_vals)=",len(all_vals))
# print("all_vals=",all_vals)
# print("all_vals.max()=",all_vals.max())
nbins = 30
nbins = 50
density_flag = False
density_flag = True
# density_flag: If False, the result will contain the number of samples in each bin. If True, the result is the value of the probability density function at the bin, normalized such that the integral over the range is 1. Note that the sum of the histogram values will not be equal to 1 unless bins of unity width are chosen; it is not a probability mass function.
# counts, bin_edges = np.histogram(all_vals, bins=nbins, density=density_flag)
         
# print("counts= ",counts)
# print("counts sum= ",counts.sum())

plt.hist(all_vals, bins=50, color='skyblue', edgecolor='black') 

# counts_pct = counts/total_count
# cdf1 = np.cumsum(counts1) * np.diff(bin_edges1)[0]
# bin_widths = np.diff(bin_edges)
# bar_centers = bin_edges[:-1] + bin_widths / 2
# plt.bar(bar_centers, counts, width=bin_widths, edgecolor='black')
# plt.bar(bar_centers, counts, width=bin_widths, edgecolor='black')
# plt.bar(bar_centers, counts_pct, width=bin_widths, edgecolor='black')

# q1 = np.percentile(all_vals, 25) # Median
q2 = np.percentile(all_vals, 50) # Median
p95 = np.percentile(all_vals, 95)

# plt.axvline(x=q2, color='b', linestyle='--', label="50th Percentile (Median)")
# plt.axvline(x=p95, color='r', linestyle='--', label="95th Percentile")


# print(all_vals[0:5])
# plt.hist(all_vals, bins=10, color='skyblue', edgecolor='black')

# plt.plot(bin_edges1[1:], cdf1, label='Sample 1 (Histogram)')

#   ------  1st plot: histogram ----------------
plt.xlabel(cell_scalar_name )
# plt.ylabel("counts (100 cells x 100 runs)")
# if output_dir_base == "cdf_1000cells":
#     plt.ylabel("counts (1000 cells x 100 runs)")
# elif output_dir_base == "cdf_200cells":
#     plt.ylabel("counts (200 cells x 100 runs)")
#plt.ylabel(f"counts ({suffix} cells x 100 runs)")
plt.ylabel(f"counts ({1000} cells x {max_runs} runs)")
# if (cell_scalar_name == "f_i"):
#     plt.xlabel("Surface fraction ($\\gamma$)")   # f_i
# elif (cell_scalar_name == "a_i"):
#     plt.xlabel("Area fraction ($\\beta$)")     # a_i

# if (cell_scalar_name == "f_i"):
#     # plt.title("PhysiCell: histogram of surface fractions ($\\beta = \\gamma = 0$)")
#     plt.title("PhysiCell: f (N=1000)")
# elif (cell_scalar_name == "a_i"):
#     # plt.title("PhysiCell: histogram of area fractions ($\\beta = \\gamma = 0$)")
#     plt.title("PhysiCell: {cell_scalar_name } (N=1000)")

# title_str = f"PhysiCell: {cell_scalar_name}  (N={max_runs})"
title_str = f"PhysiCell: {cell_scalar_name}"
plt.title(title_str)
# plt.legend()
# plt.xlim(right=1.0)
# plt.xlim(left=0.0, right=1.0)
plt.show()

#   ------  2nd plot: CDF ----------------
# plt.hist(all_vals, bins=nbins, cumulative=True, label='CDF', histtype='step', alpha=0.8, color='r')
# use "density=True"  to get normalized [0,1] y-axis
# plt.hist(all_vals, bins=100, density=True, cumulative=True, label='CDF', histtype='step', color='b')

# n, bins = np.histogram(all_vals, bins=30, density=True)
# n, bins = np.histogram(all_vals, bins=100, density=True, cumulative=True)
# plt.step(bins[:-1], n, where='pre') 
# plt.step(bins, n, where='pre') 

H,X1 = np.histogram(all_vals, bins=50, density=True)   # normed is deprecated
dx = X1[1] - X1[0]
F1 = np.cumsum(H)*dx

plt.plot(X1[1:], F1, 'b')   # smooth plot

# print("all_vals/total_count=",all_vals/total_count)
# plt.hist(all_vals/total_count, bins=nbins, cumulative=True, label='CDF', histtype='step', color='r')

# plt.hist(all_vals, bins=nbins, cumulative=True, label='CDF',  alpha=0.8, color='k')  # fill AUC
# plt.hist(all_vals, bins=nbins, cumulative=False, label='CDF', histtype='step', alpha=0.8, color='k')
# plt.axvline(x=q2, color='b', linestyle='--', label="50th Percentile (Median)")
# plt.axvline(x=p95, color='r', linestyle='--', label="95th Percentile")
#plt.title("PhysiCell: CDF of area fractions ($\\beta = \\gamma = 0$)")
# plt.ylabel("cumulative counts (1000 cells x 100 runs)")
plt.ylabel("Cumulative probability")
# if (cell_scalar_name == "f_i"):
#     # plt.title("PhysiCell: CDF of surface fractions (f_i, $\\beta = \\gamma = 0$)")
#     plt.title("PhysiCell: Free surface fraction (f_i, $\\beta = \\gamma = 0$)")
# elif (cell_scalar_name == "a_i"):
#     # plt.title("PhysiCell: CDF of area fractions ($\\beta = \\gamma = 0$)")
#     plt.title("PhysiCell: Normalized area fractions (a_i, $\\beta = \\gamma = 0$)")

plt.ylabel("CDF")
# plt.legend()
plt.xlim(right=1.0)
# plt.xlim(left=0.0, right=1.0)
# plt.xlim(left=0.0)
plt.ylim(bottom=0.0, top=1.0)
# plt.ylim(top=1.0)
# plt.ylim(bottom=0.84, top=1.0)  # for f_i
plt.grid(True)
# plt.axis('square')
# plt.axis('equal')
# fix_hist_step_vertical_line_at_end(ax)
title_str = f"PhysiCell: {cell_scalar_name}  (N={max_runs})"
plt.title(title_str)
plt.show()