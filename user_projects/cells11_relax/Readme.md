Reminder:
```
(base) M1P~/git/PhysiCell_monolayer/PhysiCell-development$ make load PROJ=cells11_relax   # or "save" if updating
(base) M1P~/git/PhysiCell_monolayer/PhysiCell-development$ make 
(base) M1P~/git/PhysiCell_monolayer/PhysiCell-development$ cp project project_11cells 

(base) M1P~/git/PhysiCell_monolayer/PhysiCell-development$ project_11cells config/monolayer_11cells_symm_repuls_10.xml 

or via Studio:
(base) M1P~/git/PhysiCell_monolayer/PhysiCell-development$ pcstudio -e project_11cells -c config/monolayer_11cells_symm_repuls_10.xml


- look for time_90pct.txt in the output folder, e.g.:
(base) M1P~/git/PhysiCell_monolayer/PhysiCell-development$ ty output_11cells_symm_repuls_10/time_90pct.txt
88.3

----------------------

- to experiment, e.g., with different repulsion values, consider using test_11cells.xml:
(base) M1P~/git/PhysiCell_monolayer/PhysiCell-development/user_projects/cells11_relax/config$ diff monolayer_11cells_symm_repuls_10.xml test_11cells.xml 
30c30
<         <folder>output_11cells_symm_repuls_10</folder>
---
>         <folder>output_test_11cells</folder>
53a54
>         <mechanics_voxel_size>40</mechanics_voxel_size>
140c141
<                     <cell_cell_repulsion_strength units="micron/min">10</cell_cell_repulsion_strength>
---
>                     <cell_cell_repulsion_strength units="micron/min">15</cell_cell_repulsion_strength>


-->
(base) M1P~/git/PhysiCell_monolayer/PhysiCell-development$ ty output_test_11cells/time_90pct.txt
58.8
----------------------

