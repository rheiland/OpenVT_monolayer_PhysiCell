
# Find the actual f_i, a_i values, given CDF percentiles

```
(base) M1P~/git/OpenVT_monolayer_PhysiCell$ python beta/chosen_CDF_pct.py 100 cdf_1000cells_linear_growth  bg00_cells1000_   f_i
-->
x values for %s: 0.00000, 0.78565, 0.89119, 0.93058, 0.94846

(base) M1P~/git/OpenVT_monolayer_PhysiCell$ python beta/chosen_CDF_pct.py 100 cdf_1000cells_linear_growth  bg00_cells1000_   a_i
-->
x values for %s: 0.49014, 0.98774, 0.99480, 0.99829, 0.99919

=== NEW: 4/23/26:
f_i:  x values for %s: 0.00000, 0.69610, 0.81921, 0.89469, 0.90317
a_i:  x values for %s: 0.29291, 0.97430, 0.98684, 0.99344, 0.99505

--> put these in param_sweep.py
f_i:  x values for %s: 0.0, 0.696, 0.819, 0.895, 0.903
a_i:  x values for %s: 0.0, 0.974, 0.987, 0.993, 0.995


```

---------------------
```
 python param_sweep.py ../project    # _no_diffusion

- when done:
 python ../beta/plot_all_new_frames.py

- while running, from another terminal:
 python ../beta/plot_cell_scalars_4states.py -s beta_or_gamma --show_colorbar -o out_cell_area_b0.9866_g0.9081 -f -1


- but still results in a lousy res .pdf
 python ../beta/plot_final_5x5_png.py
```
