'''
NOTE! be sure to use <random_seed>system_clock</random_seed> in the config file.
(base) M1P~/git/OpenVT_monolayer_PhysiCell/cdf_1000cells_linear_growth$ python ../beta/param_00_1000cells_cell_area.py ../project 


'''

To plot results in the Studio (but won't have the 4-state colormap):

'''
python ~/git/studio_dev/bin/studio.py -c bg00_cells1000_0/config.xml 
'''
and then in Plot tab, Select the output dir, then select the .mat for cells.


To plot the last frame (.mat), showing f_i or a_i::
'''
python ../beta/plot_cell_scalars-2.py -o bg00_cells1000_0 -f -1 -s a_i
or,
python ../beta/plot_cell_scalars-2.py -o bg00_cells1000_0 -f -1 -s f_i
args= Namespace(output_dir='bg00_cells1000_0', frame=-1, axes_fixed=False, colorbar_name=None, show_colorbar=False, scalar_name='f_i', xmin=None, xmax=None, ymin=None, ymax=None)
unknown= []
output_dir=bg00_cells1000_0, current_frame=-1, axes_fixed=False, colorbar=RdBu, show_colorbar=False, xmax=100.0
'''

To plot results using the 4-state colormap ("-s beta_or_gamma"; however, it's not a discrete colormap!):
'''
python ../beta/plot_cell_scalars-2.py -s beta_or_gamma --show_colorbar -o bg00_cells1000_0 -f -1 

python ../beta/plot_cell_scalars_4states.py -s beta_or_gamma --show_colorbar -o bg00_cells1000_0 -f -1

'''

<!----------------------------------------------->
To plot histograms and cumulative dist fns (CDFs), do so from the root dir:
'''

(base) M1P~/git/PhysiCell_monolayer/OpenVT_monolayer_PhysiCell$ 
python beta/all_CDF.py 100 cdf_1000cells_linear_growth bg00_cells1000_ f_i
python beta/all_CDF.py 100 cdf_1000cells_linear_growth bg00_cells1000_ a_i

python beta/all_CDF.py 1 cdf_1000cells_linear_growth bg00_cells1000_ a_i

python beta/all_CDF_percentiles.py 70 cdf_1000cells_linear_growth bg00_cells1000_ f_i

<!----------------------------------------------->
- test Dom's PDF/CDF script on a single dir:
(base) M1P~/git/PhysiCell_monolayer/OpenVT_monolayer_PhysiCell/cdf_1000cells_linear_growth$ python ../analysis/gen_csv.py bg00_cells1000_42
-->  cell_data_no_inhibition.csv
(base) M1P~/git/PhysiCell_monolayer/OpenVT_monolayer_PhysiCell/cdf_1000cells_linear_growth$ python ../analysis/dom_single.py


<!----------------------------------------------->
- prepare proper subdir and ALL .csv files for Dom's script:
(base) M1P~/git/PhysiCell_monolayer/OpenVT_monolayer_PhysiCell/cdf_1000cells_linear_growth$ 

<!----------------------------------------------->
for Dom:
python ../analysis/gen_cdf_time_series_data.py 100 5   #   <max_runs> <cell_radius>
python ../analysis/gen_cdf_idx_time_csv.py 100 443.5   #   <max_runs>  <5*T cycle duration>
mkdir PhysiCell_MonolayerGrowth_1000_Data
mv run* PhysiCell_MonolayerGrowth_1000_Data
zip -r PhysiCell_MonolayerGrowth_1000_Data.zip PhysiCell_MonolayerGrowth_1000_Data   # confirm < 100MB for github!
mv PhysiCell_MonolayerGrowth_1000_Data.zip ../results/


<!----------------------------------------------->
#---------------------------------
# the following commands now in doit.sh
python ../analysis/gen_cdf_csv.py 100   # or however many runs
mkdir PhysiCell_MonolayerGrowth_1000_Data
zip PhysiCell_MonolayerGrowth_1000_Data.zip cell_data*.csv
cp PhysiCell_MonolayerGrowth_1000_Data.zip PhysiCell_MonolayerGrowth_1000_Data
pushd PhysiCell_MonolayerGrowth_1000_Data
unzip PhysiCell_MonolayerGrowth_1000_Data.zip
mv *.zip ..
popd
python ../analysis/dom.py 100

---------------
- from this dir
python ../analysis/gen_cdf_idx_time_csv.py 100 443.5 (was 441.5, from 88.3 * 5; now using 88.7)
--> "index,t"
python ../analysis/gen_cdf_time_series_data.py 100 5.0

pushd PhysiCell_MonolayerGrowth_1000_Data
rm -rf *
mv ../run* .

Plot from:
(base) M1P~/git/doms_OpenVTMonolayerGrowthProcessing$ jn MonolayerTimeSeries_Plots.ipynb

- if it looks reasonable,
- zip up results and push to /results (if a .zip is there already, may want to save/rename it)

(base) M1P~/git/OpenVT_monolayer_PhysiCell/cdf_1000cells_linear_growth$ zip -r PhysiCell_MonolayerGrowth_1000_Data.zip PhysiCell_MonolayerGrowth_1000_Data   # all the /run dirs
(base) M1P~/git/OpenVT_monolayer_PhysiCell/cdf_1000cells_linear_growth$ mv PhysiCell_MonolayerGrowth_1000_Data.zip ../results/
(base) M1P~/git/OpenVT_monolayer_PhysiCell/cdf_1000cells_linear_growth$ pu ../results
(base) M1P~/git/OpenVT_monolayer_PhysiCell/results$ gadd PhysiCell_MonolayerGrowth_1000_Data.zip 
(base) M1P~/git/OpenVT_monolayer_PhysiCell/results$ gcommit "latest"
(base) M1P~/git/OpenVT_monolayer_PhysiCell/results$ gpush


------------------------
- Do single run of deterministic growth, i.e., standard deviation=0 for Normal Random: N(2,0)

(base) M1P~/git/OpenVT_monolayer_PhysiCell$ pcstudio -c config/deterministic_1Kcells.xml
- output_1K_deterministic
- interval = 30 min  (--> 148 frames)
- max_cells = 1024  (or just 1000)
- gamma,beta thresholds = 0,0
- cycle_duraton= 443.5
- normal_random_flag = false  (makes standard deviation=0)


(base) M1P~/git/OpenVT_monolayer_PhysiCell$ ll run_1K_deterministic|tail
...
-rw-r--r--@ 1 heiland  staff  48170 Apr 10 06:44 data000148.csv
-rw-r--r--@ 1 heiland  staff   1986 Apr 10 07:19 timeline.txt

(base) M1P~/git/OpenVT_monolayer_PhysiCell$ rm run_1K_deterministic/*

(base) M1P~/git/OpenVT_monolayer_PhysiCell$ rm output_1K_deterministic_time_series/*

# the following will append "_time_series" suffix onto the given output dir:
(base) M1P~/git/OpenVT_monolayer_PhysiCell$ python analysis/gen_time_series.py output_1K_deterministic 5
...
-->  output_1K_deterministic_time_series/data000147.csv
-->  output_1K_deterministic_time_series/data000148.csv
(base) M1P~/git/OpenVT_monolayer_PhysiCell$

(base) M1P~/git/OpenVT_monolayer_PhysiCell$ python analysis/gen_idx_time_csv.py output_1K_deterministic 443.5
sys.argv= ['analysis/gen_idx_time_csv.py', 'output_1K_deterministic', '443.5']
tau=  443.5
last_xml_file=  output_1K_deterministic/output00000148.xml
-->  output_1K_deterministic_time_series/timeline.txt

(base) M1P~/git/OpenVT_monolayer_PhysiCell$ rm -rf run_1K_deterministic
(base) M1P~/git/OpenVT_monolayer_PhysiCell$ mv output_1K_deterministic_time_series run_1K_deterministic

(base) M1P~/git/OpenVT_monolayer_PhysiCell$ zip -r run_1K_deterministic.zip run_1K_deterministic
(base) M1P~/git/OpenVT_monolayer_PhysiCell$ mv run_1K_deterministic.zip results
------------------------

'''