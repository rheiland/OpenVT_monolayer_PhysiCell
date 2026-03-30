
```
time_count.csv - single run, normalized time
actual_time_count.csv - single run, time in PhysiCell simulation mins

Metadata for the above run:
Final time for growth to 1,000 cells: 4883.0 min
dt (for force based models): 0.1
Simulation output frequency: every 50 mins (to achieve ~100 snapshots)
Prefered Area used: 78.54 (radius 5: pi * radius^2)
Growth rate used: 1 / (88.3 * 5) = 1/441.5 ~= 0.002265

In PhysiCell we compute:
    pCell->custom_data["cell_area"] =  pCell->custom_data["cell_area_0"] *
                                        (1.0 + pCell->custom_data["time_in_cycle"] / cycle_duration ); 
```
