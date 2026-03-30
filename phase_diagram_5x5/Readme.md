
# to find the actual f_i,a_i values, knowing the chosen CDF percentiles: (from root dir)
python beta/chosen_CDF_pct.py 100 cdf_1000cells_linear_growth  bg00_cells1000_   f_i    # or “a_i” as last arg

from Mar 27:
(base) M1P~/git/PhysiCell_monolayer/PhysiCell_mech_grid_xml$ python beta/chosen_CDF_pct.py 100 cdf_1000cells_linear_growth  bg00_cells1000_  f_i
-->
x values for %s: 0.00000, 0.72514, 0.83539, 0.90324, 0.91286

and:
(base) M1P~/git/PhysiCell_monolayer/PhysiCell_mech_grid_xml$ python beta/chosen_CDF_pct.py 100 cdf_1000cells_linear_growth  bg00_cells1000_  a_i
-->
x values for %s: 0.51875, 0.98036, 0.98991, 0.99496, 0.99631
---------------------

 python param_sweep.py ../project    # _no_diffusion

- when done:
 python ../beta/plot_all_new_frames.py

- while running, from another terminal:
 python ../beta/plot_cell_scalars_4states.py -s beta_or_gamma --show_colorbar -o out_cell_area_b0.9866_g0.9081 -f -1


- but still results in a lousy res .pdf
 python ../beta/plot_final_5x5_png.py 
