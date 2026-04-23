# After editing vars below (output_dir_base and output_dir) from root dir, run:
#   python beta/all_CDF.py num_runs cdf_dir dir_prefix a_i
#   e.g.,
#   python beta/all_CDF.py 100 cdf_200cells bg00_cells200_ a_i
#   python beta/all_CDF.py 100 cdf_1000cells bg00_cells1000_ a_i
#   python beta/all_CDF.py 100 cdf_1000cells_area_approx_CPM out_cells1000_ a_i
#
import sys
import os
import glob
from pyMCDS import pyMCDS
import numpy as np
import matplotlib.pyplot as plt

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
# max_runs = 10
# max_runs = 100

# output_dir_base = "cdf_100cells"
# output_dir_base = "cdf_1000cells"
# output_dir_base = "cdf_500cells_cell_area"
# output_dir_base = "cdf_1000cells_area_approx_CPM"
# output_dir_base = "cdf_200cells"

# cell_scalar_name = "f_i"
# cell_scalar_name = "a_i"
for idx in range(max_runs):  # for each of the dirs/replicate runs
    xml_file_root = "output%08d.xml" % idx
    # output_dir = output_dir_base + f'/bg00_cells100_{idx}'
    # output_dir = output_dir_base + f'/bg00_cells1000_{idx}'
    # output_dir = output_dir_base + f'/bg00_cells500_{idx}'
    # output_dir = output_dir_base + f'/out_cells1000_{idx}'
    # output_dir = output_dir_base + f'/bg00_cells200_{idx}'
    output_dir = f'{output_dir_base}/{dir_prefix}{idx}'
    print("output_dir= ",output_dir)
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
    if cell_scalar.max() < 0.4:
        print("max = ",cell_scalar.max())
    if idx == 0:
        all_vals = cell_scalar
    else:
        # all_vals = np.concatenate((all_vals, cell_scalar), axis=0)
        all_vals = np.concatenate((all_vals, cell_scalar), axis=0)
    # print(cell_scalar)

print("len(all_vals)=",len(all_vals))
# print("all_vals=",all_vals)
nbins = 10
nbins = 50
density_flag = False
# density_flag: If False, the result will contain the number of samples in each bin. If True, the result is the value of the probability density function at the bin, normalized such that the integral over the range is 1. Note that the sum of the histogram values will not be equal to 1 unless bins of unity width are chosen; it is not a probability mass function.
counts, bin_edges = np.histogram(all_vals, bins=nbins, density=density_flag)
print("counts= ",counts)
print("counts sum= ",counts.sum())
# cdf1 = np.cumsum(counts1) * np.diff(bin_edges1)[0]
bin_widths = np.diff(bin_edges)
bar_centers = bin_edges[:-1] + bin_widths / 2
plt.bar(bar_centers, counts, width=bin_widths, edgecolor='black')


# print(all_vals[0:5])
# plt.hist(all_vals, bins=10, color='skyblue', edgecolor='black')

# plt.plot(bin_edges1[1:], cdf1, label='Sample 1 (Histogram)')

plt.xlabel(cell_scalar_name )
#plt.ylabel('# of cells (avg over 100 runs)')
plt.ylabel("counts (100 cells x 100 runs)")
plt.ylabel("counts (500 cells x 100 runs)")
plt.ylabel("counts (1000 cells x 100 runs)")

if cell_scalar_name == "a_i":
    plt.xlabel("Area fraction ($\\beta$)")
    plt.title("PhysiCell: histogram of area fractions ($\\beta = \\gamma = 0$)")
elif cell_scalar_name == "f_i":
    plt.xlabel("Surface fraction ($\\gamma$)")
    plt.title("PhysiCell: histogram of surface fractions ($\\beta = \\gamma = 0$)")

plt.xlim(right=1.0)
plt.show()

#------------- 2nd plot: cumulative
if False:

    # cdf1 = np.cumsum(counts) * np.diff(bin_edges)[0]
    # plt.plot(bin_edges1[1:], cdf1, label='Sample 1 (Histogram)')
    # self.ax0.plot(bin_edges1[1:], cdf1, label='CDF')
    # plt.bar(cumulative=True, histtype='step', color='k', label='CDF'))

    #  https://matplotlib.org/stable/gallery/statistics/histogram_cumulative.html
    # plt.hist(all_vals, bins=nbins, cumulative=True, label='CDF', histtype='step', alpha=0.8, color='k')
    # plt.hist(all_vals, bins=nbins, cumulative=True, label='CDF', alpha=0.8, color='k')   # fills in solid below
    plt.ecdf(all_vals, label='CDF', alpha=0.8, color='k')   # smooth 

    # plt.xlim(right=1.0)
    plt.xlim(left=0., right=1.0)
    plt.ylim(bottom=0.5)
    plt.title("Empirical CDF")
    plt.ylabel("CDF")
    if cell_scalar_name == "a_i":
        plt.xlabel(cell_scalar_name )
    elif cell_scalar_name == "f_i":
        plt.xlabel(cell_scalar_name )
    plt.show()