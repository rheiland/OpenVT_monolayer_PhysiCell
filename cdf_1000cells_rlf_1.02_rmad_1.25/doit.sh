python ../analysis/gen_cdf_csv.py 100   # or however many runs
mkdir PhysiCell_MonolayerGrowth_1000_Data
zip PhysiCell_MonolayerGrowth_1000_Data.zip cell_data*.csv
cp PhysiCell_MonolayerGrowth_1000_Data.zip PhysiCell_MonolayerGrowth_1000_Data
pushd PhysiCell_MonolayerGrowth_1000_Data
unzip PhysiCell_MonolayerGrowth_1000_Data.zip
mv *.zip ..
popd
python ../analysis/dom.py 100
