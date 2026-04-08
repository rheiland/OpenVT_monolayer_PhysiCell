import sys
import glob
import argparse
import string
import os
import time
import xml.etree.ElementTree as ET  # https://docs.python.org/2/library/xml.etree.elementtree.html
from pathlib import Path

import matplotlib
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle, Ellipse, Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.colors as mplc
import cmaps
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib import gridspec
from collections import deque
import math
import numpy as np
# import scipy.io  # .io.loadmat(full_fname, info_dict)
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from pyMCDS import pyMCDS

# argc = len(sys.argv)
# if argc < 5:
#     print("Usage:  args = output_dir colorbar_name frame show_colorbar [xmin xmax ymin ymax]")
#     sys.exit()
# print('argv=',sys.argv)
# print('argc = len(argv)=',len(sys.argv))

# idx = 0
# print('argv[0]=',sys.argv[idx])

# idx += 1
# output_dir = sys.argv[idx]
# idx += 1
# colorbar_name = sys.argv[idx]
# print("colorbar_name= ",colorbar_name)
# idx += 1
# current_frame = int(sys.argv[idx])
# print("current_frame= ",current_frame)
# # p1=string.atof(sys.argv[1])

# idx += 1
# show_colorbar = int(sys.argv[idx])

# print("idx, argc= ",idx,argc)
# if idx + 1 == argc:
#     fixed_axes = False
# else:
#     fixed_axes = True
#     idx += 1
#     plot_xmin = float(sys.argv[idx])
#     idx += 1
#     plot_xmax = float(sys.argv[idx])
#     idx += 1
#     plot_ymin = float(sys.argv[idx])
#     idx += 1
#     plot_ymax = float(sys.argv[idx])


class Vis():
    def __init__(self, output_dir, current_frame,axes_fixed,colorbar_name,show_colorbar,scalar_name,plot_xmin,plot_xmax,plot_ymin,plot_ymax):
        super().__init__()
        # global self.config_params

        print("\n--- Vis()  plot_xmax=", plot_xmax)

        self.xml_root = None
        self.current_frame = current_frame
        # self.timer = QtCore.QTimer()
        # self.t.timeout.connect(self.task)
        # self.timer.timeout.connect(self.play_plot_cb)

        # self.tab = QWidget()
        # self.tabs.resize(200,5)
        
        self.use_defaults = True
        self.title_str = ""
        self.xmin = -1000
        self.xmax = 1000
        self.x_range = self.xmax - self.xmin

        self.ymin = -1000
        self.ymax = 1000
        self.y_range = self.ymax - self.ymin
        self.show_nucleus = False
        self.show_edge = False
        self.alpha = 0.7
        # self.cells_toggle = None
        # self.substrates_toggle = None

        basic_length = 12.0
        # rwh: no longer relevant??
        self.figsize_width_substrate = 15.0  # allow extra for colormap
        self.figsize_height_substrate = basic_length

        # self.figsize_width_2Dplot = basic_length
        # self.figsize_height_2Dplot = basic_length

        # rwh: TODO these params
        self.modulo = 1
        self.field_index = 4
        # define dummy size of mesh (set in the tool's primary module)
        self.numx = 0
        self.numy = 0
        self.colormap_min = 0.0
        self.colormap_max = 10.0
        self.colormap_fixed_toggle = False
        self.fontsize = 10


        # self.width_substrate = basic_length  # allow extra for colormap
        # self.height_substrate = basic_length

        self.figsize_width_svg = basic_length
        self.figsize_height_svg = basic_length

        # self.output_dir = "/Users/heiland/dev/PhysiCell_V.1.8.0_release/output"
        self.output_dir = "./output"
        self.output_dir = output_dir

        self.fixed_axes = axes_fixed
        if self.fixed_axes:
            self.plot_xmin = plot_xmin
            self.plot_xmax = plot_xmax
            self.plot_ymin = plot_ymin
            self.plot_ymax = plot_ymax
            print("\n--- setting  plot_xmax=", self.plot_xmax)

        self.cbar_name = colorbar_name

        self.show_colorbar = show_colorbar
        self.scalar_name = scalar_name

        self.customized_output_freq = False

        # rwh hacks (originally gen'd vis_base.py)
        self.cell_scalar_updated = True
        self.view_aspect_square = True
        self.cbar_label_fontsize = 8
        self.cbar_label_fontsize = 10
        self.title_fontsize = 10
        self.cbar2 = None
        self.figure = None
        self.cax2 = None
        self.fix_cells_cmap_flag = False
        self.cell_fill = True
        # self.cell_edge = True
        self.cell_edge = False
        self.cell_line_width = 0.5
        self.discrete_cell_scalars = []
        self. cell_scalar_human2mcds_dict=  {'a_i': 'a_i', 'apoptotic_phagocytosis_rate': 'apoptotic_phagocytosis_rate', 'asymmetric_division_probabilities': 'asymmetric_division_probabilities', 'attachment_elastic_constant': 'attachment_elastic_constant', 'attachment_rate': 'attachment_rate', 'attack_damage_rate': 'attack_damage_rate', 'attack_duration': 'attack_duration', 'attack_rates': 'attack_rates', 'attack_target': 'attack_target', 'attack_total_damage_delivered': 'attack_total_damage_delivered', 'beta_or_gamma': 'beta_or_gamma', 'calcification_rate': 'calcification_rate', 'calcified_fraction': 'calcified_fraction', 'cell_BM_adhesion_strength': 'cell_BM_adhesion_strength', 'cell_BM_repulsion_strength': 'cell_BM_repulsion_strength', 'cell_adhesion_affinities': 'cell_adhesion_affinities', 'cell_cell_adhesion_strength': 'cell_cell_adhesion_strength', 'cell_cell_repulsion_strength': 'cell_cell_repulsion_strength', 'cell_type': 'cell_type', 'chemotactic_sensitivities': 'chemotactic_sensitivities', 'chemotaxis_direction': 'chemotaxis_direction', 'chemotaxis_index': 'chemotaxis_index', 'contact_with_basement_membrane': 'contact_with_basement_membrane', 'current_cycle_phase_exit_rate': 'current_cycle_phase_exit_rate', 'current_death_model': 'current_death_model', 'current_phase': 'current_phase', 'cycle_model': 'cycle_model', 'cytoplasmic_biomass_change_rate': 'cytoplasmic_biomass_change_rate', 'cytoplasmic_volume': 'cytoplasmic_volume', 'damage': 'damage', 'damage_rate': 'damage_rate', 'damage_repair_rate': 'damage_repair_rate', 'dead': 'dead', 'death_rates_0': 'death_rates_0', 'death_rates_1': 'death_rates_1', 'detachment_rate': 'detachment_rate', 'elapsed_time_in_phase': 'elapsed_time_in_phase', 'f_i': 'f_i', 'fluid_change_rate': 'fluid_change_rate', 'fluid_fraction': 'fluid_fraction', 'fraction_released_at_death': 'fraction_released_at_death', 'fraction_transferred_when_ingested': 'fraction_transferred_when_ingested', 'fusion_rates': 'fusion_rates', 'immunogenicities': 'immunogenicities', 'internalized_total_substrates': 'internalized_total_substrates', 'is_motile': 'is_motile', 'live_phagocytosis_rates': 'live_phagocytosis_rates', 'maximum_number_of_attachments': 'maximum_number_of_attachments', 'migration_bias': 'migration_bias', 'migration_bias_direction_x': 'migration_bias_direction_x', 'migration_bias_direction_y': 'migration_bias_direction_y', 'migration_bias_direction_z': 'migration_bias_direction_z', 'migration_speed': 'migration_speed', 'motility_vector_x': 'motility_vector_x', 'motility_vector_y': 'motility_vector_y', 'motility_vector_z': 'motility_vector_z', 'necrotic_phagocytosis_rate': 'necrotic_phagocytosis_rate', 'net_export_rates': 'net_export_rates', 'nuclear_biomass_change_rate': 'nuclear_biomass_change_rate', 'nuclear_radius': 'nuclear_radius', 'nuclear_volume': 'nuclear_volume', 'number_of_nuclei': 'number_of_nuclei', 'orientation_x': 'orientation_x', 'orientation_y': 'orientation_y', 'orientation_z': 'orientation_z', 'other_dead_phagocytosis_rate': 'other_dead_phagocytosis_rate', 'persistence_time': 'persistence_time', 'polarity': 'polarity', 'position_x': 'position_x', 'position_y': 'position_y', 'position_z': 'position_z', 'pressure': 'pressure', 'radius': 'radius', 'relative_maximum_adhesion_distance': 'relative_maximum_adhesion_distance', 'saturation_densities': 'saturation_densities', 'secretion_rates': 'secretion_rates', 'surface_area': 'surface_area', 'target_fluid_fraction': 'target_fluid_fraction', 'target_solid_cytoplasmic': 'target_solid_cytoplasmic', 'target_solid_nuclear': 'target_solid_nuclear', 'total_attack_time': 'total_attack_time', 'total_volume': 'total_volume', 'transformation_rates': 'transformation_rates', 'uptake_rates': 'uptake_rates', 'velocity_x': 'velocity_x', 'velocity_y': 'velocity_y', 'velocity_z': 'velocity_z', '(probability of) asymmetric division to default': 'asymmetric_division_probabilities', '(rate of) attack default': 'attack_rates', 'adhesive affinity to default': 'cell_adhesion_affinities', 'chemotactic response to substrate': 'chemotactic_sensitivities', 'fraction released at death of substrate': 'fraction_released_at_death', 'fraction transferred when ingested of substrate': 'fraction_transferred_when_ingested', '(rate of) fuse to default': 'fusion_rates', 'immunogenicity to default': 'immunogenicities', '(amount of) intracellular substrate': 'internalized_total_substrates', '(rate of) phagocytose default': 'live_phagocytosis_rates', '(rate of) substrate export': 'net_export_rates', 'substrate secretion target': 'saturation_densities', '(rate of) substrate secretion ': 'secretion_rates', '(rate of) transform to default': 'transformation_rates', '(rate of) substrate uptake': 'uptake_rates'}

        #-------------------------------------------
        label_width = 110
        domain_value_width = 100
        value_width = 60
        label_height = 20
        units_width = 70

        self.create_figure()

        self.reset_plot_cb("")

        if self.current_frame < 0:
            xml_pattern = self.output_dir + "/" + "output*.xml"
            xml_files = glob.glob(xml_pattern)
            xml_files.sort()
            # print("xml_files= ",xml_files)
            last_file = xml_files[-1]
            print("last file= ",last_file)
            print("last index= ",last_file[-12:-4])
            self.current_frame = int(last_file[-12:-4])
            print("current_frame = ",self.current_frame)
            # self.current_frame = 3

        # try:
            self.plot_cell_scalar(self.current_frame)
            png_filename = Path(self.output_dir,f'frame{self.current_frame}')
            plt.savefig(png_filename)
            png_filename = Path(self.output_dir,'keep.png')
            # plt.savefig(png_filename)
            # plt.savefig(png_filename, bbox_inches='tight', pad_inches=0)
            plt.savefig(png_filename, bbox_inches='tight')

            # plt.show()   # rwh - comment out if just generating 5x5 keep.png
        # except:
            # print("-- error plotting")


    def get_mcds_cells_df(self, mcds):
        try:
            df_all_cells = mcds.get_cell_df()
        except:
            print("vis_tab.py: plot_cell_scalar(): error performing mcds.get_cell_df()")
            return
        # if self.celltype_filter:  #rwh
        if False:
            return df_all_cells.loc[ df_all_cells['cell_type'].isin(self.celltype_filter) ]
        else:
            return df_all_cells
        

    def reset_plot_cb(self, text):
        # print("-------------- reset_plot_cb() ----------------")
        xml_file = Path(self.output_dir, "initial.xml")
        if not os.path.isfile(xml_file):
            print("Expecting initial.xml, but does not exist.")
            return


        tree = ET.parse(Path(self.output_dir, "initial.xml"))
        xml_root = tree.getroot()

        bds_str = xml_root.find(".//microenvironment//domain//mesh//bounding_box").text
        bds = bds_str.split()
        print('bds=',bds)
        self.xmin = float(bds[0])
        self.xmax = float(bds[3])
        self.x_range = self.xmax - self.xmin

        self.ymin = float(bds[1])
        self.ymax = float(bds[4])
        self.y_range = self.ymax - self.ymin

        # self.numx =  math.ceil( (self.xmax - self.xmin) / config_tab.xdelta.value)
        # self.numy =  math.ceil( (self.ymax - self.ymin) / config_tab.ydelta.value)
        self.numx =  math.ceil( (self.xmax - self.xmin) / 20.)
        self.numy =  math.ceil( (self.ymax - self.ymin) / 20.)
        # print(" calc: numx,numy = ",self.numx, self.numy)

        # self.current_frame = current_frame
        print('frame # ',self.current_frame)
        self.plot_cell_scalar(self.current_frame)
        # self.plot_substrate(self.current_frame)
        # self.canvas.update()
        # self.canvas.draw()

    #---------------------------------------------------------------------------
    def create_figure(self):
        # print("\n--   vis_tab_cell_scalars.py: --------- create_figure(): ------- creating figure, canvas, ax0")
        if self.figure is not None:
            print("              self.figure is None, so return!")
            return
        # self.figure = plt.figure()
        self.figure = plt.figure(figsize=(5, 5.5))
        # self.gs = gridspec.GridSpec(2,2, height_ratios=[20,1], width_ratios=[20,1]) # top row is [plot, substrate colorbar]; bottom row is [cells colorbar, nothing]
        self.gs = gridspec.GridSpec(2,2, height_ratios=[20,0.8], width_ratios=[20,0.5]) # top row is [plot, substrate colorbar]; bottom row is [cells colorbar, nothing]
        # self.canvas = FigureCanvasQTAgg(self.figure)
        # print("     self.canvas= ",self.canvas)
        # self.canvas.setStyleSheet("background-color:transparent;")

        # Adding one subplot for image
        # self.ax0 = self.figure.add_subplot(111)
        # self.ax0 = self.figure.add_subplot(111, adjustable='box', aspect=1.2)
        # self.ax0 = self.figure.add_subplot(111, adjustable='box', aspect=self.aspect_ratio)
        self.ax0 = self.figure.add_subplot(self.gs[0,0], adjustable='box')
        
        # self.ax0.get_xaxis().set_visible(False)
        # self.ax0.get_yaxis().set_visible(False)
        # plt.tight_layout()

        # self.reset_model()  # rwh

    #---------------------------------------------------------------------------
    def circles(self, x, y, s, c='b', vmin=None, vmax=None, **kwargs):
        """
        See https://gist.github.com/syrte/592a062c562cd2a98a83 

        Make a scatter plot of circles. 
        Similar to plt.scatter, but the size of circles are in data scale.
        Parameters
        ----------
        x, y : scalar or array_like, shape (n, )
            Input data
        s : scalar or array_like, shape (n, ) 
            Radius of circles.
        c : color or sequence of color, optional, default : 'b'
            `c` can be a single color format string, or a sequence of color
            specifications of length `N`, or a sequence of `N` numbers to be
            mapped to colors using the `cmap` and `norm` specified via kwargs.
            Note that `c` should not be a single numeric RGB or RGBA sequence 
            because that is indistinguishable from an array of values
            to be colormapped. (If you insist, use `color` instead.)  
            `c` can be a 2-D array in which the rows are RGB or RGBA, however. 
        vmin, vmax : scalar, optional, default: None
            `vmin` and `vmax` are used in conjunction with `norm` to normalize
            luminance data.  If either are `None`, the min and max of the
            color array is used.
        kwargs : `~matplotlib.collections.Collection` properties
            Eg. alpha, edgecolor(ec), facecolor(fc), linewidth(lw), linestyle(ls), 
            norm, cmap, transform, etc.
        Returns
        -------
        paths : `~matplotlib.collections.PathCollection`
        Examples
        --------
        a = np.arange(11)
        circles(a, a, s=a*0.2, c=a, alpha=0.5, ec='none')
        plt.colorbar()
        License
        --------
        This code is under [The BSD 3-Clause License]
        (http://opensource.org/licenses/BSD-3-Clause)
        """

        if np.isscalar(c):
            kwargs.setdefault('color', c)
            c = None

        if 'fc' in kwargs:
            kwargs.setdefault('facecolor', kwargs.pop('fc'))
        if 'ec' in kwargs:
            kwargs.setdefault('edgecolor', kwargs.pop('ec'))
        if 'ls' in kwargs:
            kwargs.setdefault('linestyle', kwargs.pop('ls'))
        if 'lw' in kwargs:
            kwargs.setdefault('linewidth', kwargs.pop('lw'))
        # You can set `facecolor` with an array for each patch,
        # while you can only set `facecolors` with a value for all.

        zipped = np.broadcast(x, y, s)
        patches = [Circle((x_, y_), s_)
                for x_, y_, s_ in zipped]
        collection = PatchCollection(patches, **kwargs)
        if c is not None:
            # print("--- circles(): type(c)=",type(c))
            c = c.values
            # print("--- circles() (2): type(c)=",type(c))

            # print("--- circles(): c=",c)
            c = np.broadcast_to(c, zipped.shape).ravel()
            collection.set_array(c)
            # print("--- circles(): vmin,vmax=",vmin,vmax)
            collection.set_clim(vmin, vmax)

        # ax = plt.gca()
        # ax.add_collection(collection)
        # ax.autoscale_view()
        self.ax0.add_collection(collection)
        self.ax0.autoscale_view()
        plt.draw_if_interactive()
        # if c is not None:
        #     try:
        #         print("------ circles(): doing plt.sci(collection), type(collection)=",type(collection))
        #         plt.sci(collection)
        #         # self.ax0.sci(collection)
        #         # self.ax0.sci(collection)
        #     except:
        #         print("--- ERROR in circles() doing plt.sci(collection)")
        return collection

    #------------------------------------------------------------
    def plot_cell_scalar(self, frame):
        print("~~~ plot_cell_scalar()")
        self.disable_cell_scalar_cb = False
        if self.disable_cell_scalar_cb:
            return
            
        # if self.show_voxel_grid:
        #     self.plot_voxel_grid()
        # if self.show_mechanics_grid:
        #     self.plot_mechanics_grid()


        xml_file_root = "output%08d.xml" % frame
        xml_file_root = "output%08d.xml" % self.current_frame
        xml_file = os.path.join(self.output_dir, xml_file_root)
        # cell_scalar_humanreadable_name = self.cell_scalar_combobox.currentText()
        # cell_scalar_humanreadable_name = "a_i"
        cell_scalar_humanreadable_name = self.scalar_name
        if cell_scalar_humanreadable_name in self.cell_scalar_human2mcds_dict.keys():
            cell_scalar_mcds_name = self.cell_scalar_human2mcds_dict[cell_scalar_humanreadable_name]
        else:
            cell_scalar_mcds_name = cell_scalar_humanreadable_name
        # cbar_name = self.cell_scalar_cbar_combobox.currentText()
        # cbar_name = "RdBu"
        cbar_name = self.cbar_name
        if not Path(xml_file).is_file():
            print("ERROR: file not found",xml_file)
            return

        mcds = pyMCDS(xml_file_root, self.output_dir, microenv=False, graph=False, verbose=False)
        df_cells = self.get_mcds_cells_df(mcds)
        total_min = mcds.get_time()  # warning: can return float that's epsilon from integer value

        try:
            cell_scalar = df_cells[cell_scalar_mcds_name]
        except:
            print("vis_tab.py: plot_cell_scalar(): error performing df_cells[cell_scalar_mcds_name]")

            # self.cells_svg_rb.setChecked(True)
            # self.plot_cells_svg = True
            # self.disable_cell_scalar_widgets()
            return
                
        if self.fix_cells_cmap_flag:
            vmin = self.cells_cmin_value
            vmax = self.cells_cmax_value
        else:
            vmin = cell_scalar.min()
            vmax = cell_scalar.max()
        # print("---- vmin,vmax= ",vmin,vmax)

        num_cells = len(cell_scalar)
        # print("  len(cell_scalar) = ",len(cell_scalar))
        # fix_cmap = 0
        # print(f'   cell_scalar.min(), max() = {vmin}, {vmax}')
        cell_vol = df_cells['total_volume']
        # print(f'   cell_vol.min(), max() = {cell_vol.min()}, {cell_vol.max()}')

        four_thirds_pi =  4.188790204786391
        cell_radii = np.divide(cell_vol, four_thirds_pi)
        cell_radii = np.power(cell_radii, 0.333333333333333333333333333333333333333)

        xvals = df_cells['position_x']
        yvals = df_cells['position_y']

        # self.title_str += "   cells: " + svals[2] + "d, " + svals[4] + "h, " + svals[7][:-3] + "m"
        # self.title_str = "(" + str(frame) + ") Current time: " + str(total_min) + "m"
        
        #-----------------------------------------------------
        # rwh - unique to OpenVT monolayer model
        if cell_scalar_mcds_name == "beta_or_gamma": 
            print("------- doing discrete beta_or_gamma")
            self.discrete_variable = [0,1,2,3]

            # names_observed = ["phase #%d" % i for i in sorted(list(self.discrete_variable_observed)) if i in [0,1,2,3]]
            from_list = matplotlib.colors.LinearSegmentedColormap.from_list
            self.discrete_variable.sort()
            if (len(self.discrete_variable) == 1): 
                cbar_name = from_list(None, cmaps.gray_gray[0:2], len(self.discrete_variable))  # annoying hack
            else: 
                try:
                    # cbar_name = from_list(None, cmaps.paint_clist[0:len(self.discrete_variable)], len(self.discrete_variable))
                    # Lutz: light-green(uninhibited), light-blue, yellow, red
                    # cbar_name = from_list(None, [[0.5, 1, 0.5],[0,0.5,1],[1,1,0],[1,0,0]], len(self.discrete_variable))
                    # Roman's cartography app: dark blue, light-blue, yellow, red
                    # Red: 215,25,28
                    # yellow: 253,174,97
                    # light blue: 171,217,233
                    # blue: 44,123,182

                    # cbar_name = from_list(None, [[44./255, 123./255, 182./255],[171./255,217./255,233./255],[253./255,174./255,97./255],[215./255,25./255,28./255]], len(self.discrete_variable))
                    cbar_name = from_list(None, [[44/255.,123/255.,182/255.], [253/255.,174/255.,97/255.], 
                        [171/255.,217/255.,233/255.],[215/255.,25/255.,28/255.]], len(self.discrete_variable))
                    # print("cmaps.paint_clist=",cmaps.paint_clist)
                    # print("cbar_name=",cbar_name) # <matplotlib.colors.LinearSegmentedColormap
                except:
                    return

            # usual categorical colormap on matplotlib has at max 20 colors (using colorcet the colormap glasbey_bw has n colors )
            # cbar_name = from_list(None, cc.glasbey_bw, len(self.discrete_variable))
            vmin = None
            vmax = None
            # Change the values between 0 and number of possible values
            for i, value in enumerate(self.discrete_variable):
                cell_scalar = cell_scalar.replace(value,i)

        #-----------------------------------------------------
        if cell_scalar_mcds_name in self.discrete_cell_scalars: 

            self.discrete_variable_observed = self.discrete_variable_observed.union(set([int(i) for i in np.unique(cell_scalar)]))

            if cell_scalar_mcds_name == "current_phase":   # and if "Fixed" range is checked
                self.discrete_variable = list(self.cycle_phases.keys())
                names_observed = [self.cycle_phases[i] for i in sorted(list(self.discrete_variable_observed)) if i in self.cycle_phases.keys()]

            elif cell_scalar_mcds_name == "cell_type":
                # I'm not sure I should be calling this every time. But I'm also not sure about the life cycle of celltype_name
                self.get_cell_types_from_config()
                self.discrete_variable = list(range(len(self.celltype_name)))
                names_observed = [self.celltype_name[i] for i in sorted(list(self.discrete_variable_observed)) if i < len(self.celltype_name)]
                
            elif cell_scalar_mcds_name == "cycle_model":
                self.discrete_variable = list(self.cycle_models.keys())
                names_observed = [self.cycle_models[i] for i in sorted(list(self.discrete_variable_observed)) if i in self.cycle_models.keys()]

            elif cell_scalar_mcds_name == "current_death_model":
                self.discrete_variable = [0,1]
                names_observed = ["phase #%d" % i for i in sorted(list(self.discrete_variable_observed)) if i in [0,1]]
            
            elif cell_scalar_mcds_name == "is_motile":
                self.discrete_variable = [0,1]
                names_observed = ["motile" if i == 1 else "stationnary" for i in sorted(list(self.discrete_variable_observed)) if i in [0,1]]
                
            elif cell_scalar_mcds_name == "dead":
                self.discrete_variable = [0,1]
                names_observed = ["dead" if i == 1 else "alive" for i in sorted(list(self.discrete_variable_observed)) if i in [0,1]]


            else:
                self.discrete_variable = [int(i) for i in list(set(cell_scalar))] # It's a set of possible value of the variable
                names_observed = [str(int(i)) for i in sorted(list(self.discrete_variable_observed))] 

        # if( discrete_variable ): # Generic way: if variable is discrete
            # self.cell_scalar_cbar_combobox.setEnabled(False)
            from_list = matplotlib.colors.LinearSegmentedColormap.from_list
            self.discrete_variable.sort()
            if (len(self.discrete_variable) == 1): 
                cbar_name = from_list(None, cmaps.gray_gray[0:2], len(self.discrete_variable))  # annoying hack
            else: 
                try:
                    cbar_name = from_list(None, cmaps.paint_clist[0:len(self.discrete_variable)], len(self.discrete_variable))
                except:
                    return

            # usual categorical colormap on matplotlib has at max 20 colors (using colorcet the colormap glasbey_bw has n colors )
            # cbar_name = from_list(None, cc.glasbey_bw, len(self.discrete_variable))
            vmin = None
            vmax = None
            # Change the values between 0 and number of possible values
            for i, value in enumerate(self.discrete_variable):
                cell_scalar = cell_scalar.replace(value,i)
                # print("cell_scalar=",cell_scalar)
        else: 
            # self.cell_scalar_cbar_combobox.setEnabled(True)
            self.discrete_variable = None
            self.discrete_variable_observed = set()
            
        mins = round(total_min)  # hack, assume we want integer mins
        hrs = int(mins/60)
        days = int(hrs/24)
        # print(f"mins={mins}, hrs={hrs}, days={days}")
        self.title_str = '%d days, %d hrs, %d mins' % (days, hrs-days*24, mins-hrs*60)
        # self.title_str = '%f mins' % (total_min)  # rwh: custom
        # self.title_str += " (" + str(num_cells) + " agents)"
        self.title_str += " (" + str(num_cells) + " cells)"

        # self.title_str = f'{self.output_dir} 10K, {mins} mins'   # rwh: for phase diagram 5x5

        axes_min = mcds.get_mesh()[0][0][0][0]
        axes_max = mcds.get_mesh()[0][0][-1][0]

        if (self.cell_fill):
            if (self.cell_edge):
                try:
                    # print("plot circles with vmin,vmax=",vmin,vmax)  # None,None
                    # print("plot circles with cbar_name=",cbar_name)  # <matplotlib.colors.LinearSegmentedColormap object at 0x1690d5330>
                    cell_plot = self.circles(xvals,yvals, s=cell_radii, c=cell_scalar, edgecolor='black', linewidth=self.cell_line_width, cmap=cbar_name, vmin=vmin, vmax=vmax)
                    # cell_plot = self.circles(xvals,yvals, s=cell_radii, edgecolor=cell_scalar, linewidth=0.5, cmap=cbar_name, vmin=vmin, vmax=vmax)
                except (ValueError):
                    print("\n------ ERROR: Exception from circles with edges\n")
                    pass
            else:
                # cell_plot = self.circles(xvals,yvals, s=cell_radii, c=cell_scalar, cmap=cbar_name)
                cell_plot = self.circles(xvals,yvals, s=cell_radii, c=cell_scalar, cmap=cbar_name, vmin=vmin, vmax=vmax)

        else:  # semi-trransparent cell, but with (thicker) edge  (TODO: how to make totally transparent?)
            if (self.cell_edge):
                cell_plot = self.circles(xvals,yvals, s=cell_radii, c=cell_scalar, edgecolor='black', linewidth=self.cell_line_width2, cmap=cbar_name, vmin=vmin, vmax=vmax, alpha=self.cell_alpha)
            else:
                cell_plot = self.circles(xvals,yvals, s=cell_radii, c=cell_scalar, cmap=cbar_name, vmin=vmin, vmax=vmax, alpha=self.cell_alpha)


        # print("------- plot_cell_scalar() -------------")
        num_axes =  len(self.figure.axes)
        # print("# axes = ",num_axes)
        # if num_axes > 1: 
        # if self.axis_id_cellscalar:
    
        if self.show_colorbar:

            # rwh:
            if( self.discrete_variable ): # Generic way: if variable is discrete
            # if( self.discrete_variable or cell_scalar_mcds_name == "beta_or_gamma"): # Generic way: if variable is discrete
                # Then we don't need the cax2
                if self.cax2 is not None:
                    try:
                        self.cax2.remove()
                        self.cax2 = None
                    except:
                        pass
                # Coloring the cells as it used to be
                cell_plot.set_clim(vmin=-0.5,vmax=len(self.discrete_variable)-0.5) 
                
                # Creating empty plots to add the legend
                lp = lambda i: plt.plot([],color=cmaps.paint_clist[i], ms=np.sqrt(81), mec="none",
                                        label="Feature {:g}".format(i), ls="", marker="o")[0]
                handles = [lp(self.discrete_variable.index(i)) for i in sorted(list(self.discrete_variable_observed)) if i in self.discrete_variable]
                try: # cautionary for out of date mpl versions, e.g., nanoHUB
                    self.ax0.legend(handles=handles,title=cell_scalar_humanreadable_name, labels=names_observed, loc='upper center', bbox_to_anchor=(0.5, -0.15),ncols=4)
                except:
                    pass

            else:   # Note: vis_tab_ecm.py seems to avoid any memory leak and with simpler code
                # If it's not there, we create it
                if self.cax2 is None:
                    self.cax2 = self.figure.add_subplot(self.gs[1,0])
                    # ax2_divider = make_axes_locatable(self.ax0)
                    # self.cax2 = ax2_divider.append_axes("bottom", size="4%", pad="8%")
                if self.cbar2 is None:
                    self.cbar2 = self.figure.colorbar(cell_plot, ticks=None, cax=self.cax2, orientation="horizontal")
                    self.cbar2.ax.tick_params(labelsize=self.fontsize)
                elif self.cell_scalar_updated:
                    print("cbar #2")
                    self.cbar2 = self.figure.colorbar(cell_plot, ticks=None, cax=self.cax2, orientation="horizontal")
                    self.cell_scalar_updated = False
                else:
                    self.cbar2.update_normal(cell_plot)  # partial fix for memory leak

                # self.cbar_label_fontsize = 0  #rwh
                if cell_scalar_mcds_name != "beta_or_gamma": 
                    self.cbar2.ax.set_xlabel(cell_scalar_humanreadable_name, fontsize=self.cbar_label_fontsize)
   
        self.ax0.set_title(self.title_str, fontsize=self.title_fontsize)
        self.ax0.set_xlabel(self.output_dir)

        # rwh
        if self.fixed_axes:
            print("plotting fixed_axes:  xmax=", self.plot_xmax)
            self.ax0.set_xlim(self.plot_xmin, self.plot_xmax)
            self.ax0.set_ylim(self.plot_ymin, self.plot_ymax)

        if self.view_aspect_square:
            self.ax0.set_aspect('equal')
        else:
            self.ax0.set_aspect('auto')

    #---------------------------------------------------------------------------
    # assume "frame" is cell frame #, unless Cells is togggled off, then it's the substrate frame #
    # def plot_substrate(self, frame, grid):
    def plot_substrate_old(self, frame):   # rwh: no longer used

        # print("plot_substrate(): frame*self.substrate_delta_t  = ",frame*self.substrate_delta_t)
        # print("plot_substrate(): frame*self.svg_delta_t  = ",frame*self.svg_delta_t)
        # print("plot_substrate(): fig width: SVG+2D = ",self.figsize_width_svg + self.figsize_width_2Dplot)  # 24
        # print("plot_substrate(): fig width: substrate+2D = ",self.figsize_width_substrate + self.figsize_width_2Dplot)  # 27

        self.title_str = ''

        # Recall:
        # self.svg_delta_t = config_tab.svg_interval.value
        # self.substrate_delta_t = config_tab.mcds_interval.value
        # self.modulo = int(self.substrate_delta_t / self.svg_delta_t)
        # self.therapy_activation_time = user_params_tab.therapy_activation_time.value

        # print("plot_substrate(): pre_therapy: max svg, substrate frames = ",max_svg_frame_pre_therapy, max_substrate_frame_pre_therapy)

        # Assume: # .svg files >= # substrate files
#        if (self.cells_toggle.value):

        if self.substrates_toggle.isChecked():
            self.fig, (self.ax0) = plt.subplots(1, 1, figsize=(self.figsize_width_substrate, self.figsize_height_substrate))

            self.substrate_frame = int(frame / self.modulo)

            fname = "output%08d_microenvironment0.mat" % self.substrate_frame
            xml_fname = "output%08d.xml" % self.substrate_frame
            # fullname = output_dir_str + fname

    #        fullname = fname
            full_fname = os.path.join(self.output_dir, fname)
            print("--- plot_substrate(): full_fname=",full_fname)
            full_xml_fname = os.path.join(self.output_dir, xml_fname)
    #        self.output_dir = '.'

    #        if not os.path.isfile(fullname):
            if not os.path.isfile(full_fname):
                print("Once output files are generated, click the slider.")  # No:  output00000000_microenvironment0.mat
                return

    #        tree = ET.parse(xml_fname)
            tree = ET.parse(full_xml_fname)
            xml_root = tree.getroot()
            mins = round(int(float(xml_root.find(".//current_time").text)))  # TODO: check units = mins
            self.substrate_mins= round(int(float(xml_root.find(".//current_time").text)))  # TODO: check units = mins

            hrs = int(mins/60)
            days = int(hrs/24)
            self.title_str = 'substrate: %dd, %dh, %dm' % (int(days),(hrs%24), mins - (hrs*60))
            # self.title_str = 'substrate: %dm' % (mins )   # rwh

            info_dict = {}
            scipy.io.loadmat(full_fname, info_dict)
            M = info_dict['multiscale_microenvironment']
            f = M[self.field_index, :]   # 4=tumor cells field, 5=blood vessel density, 6=growth substrate

            try:
                print("numx, numy = ",self.numx, self.numy)
                xgrid = M[0, :].reshape(self.numy, self.numx)
                ygrid = M[1, :].reshape(self.numy, self.numx)
            except:
                print("substrates.py: mismatched mesh size for reshape: numx,numy=",self.numx, self.numy)
                pass
#                xgrid = M[0, :].reshape(self.numy, self.numx)
#                ygrid = M[1, :].reshape(self.numy, self.numx)

            num_contours = 15
            # levels = MaxNLocator(nbins=num_contours).tick_values(self.colormap_min.value, self.colormap_max.value)
            levels = MaxNLocator(nbins=num_contours).tick_values(self.colormap_min, self.colormap_max)
            contour_ok = True
            # if (self.colormap_fixed_toggle.isChecked()):
            if (self.colormap_fixed_toggle):
                try:
                    substrate_plot = self.ax0.contourf(xgrid, ygrid, M[self.field_index, :].reshape(self.numy, self.numx), 
                            levels=levels, extend='both', cmap="viridis", fontsize=self.fontsize)
                except:
                    contour_ok = False
                    print('got error on contourf 1.')
            else:    
                try:
                    print("field min,max= ", M[self.field_index, :].min(), M[self.field_index, :].max())
                    substrate_plot = self.ax0.contourf(xgrid, ygrid, M[self.field_index, :].reshape(self.numy,self.numx), 
                            num_contours, cmap = "viridis" ) #  cmap=self.colormap_dd.value)
                except:
                    contour_ok = False
                    print('\n -->> got error on contourf 2  \n')  # rwh: argh, getting here

            if (contour_ok):
                self.ax0.set_title(self.title_str, fontsize=self.fontsize)
                # cbar = self.figure.colorbar(substrate_plot, ax=self.ax0)
                cbar = self.figure.colorbar(substrate_plot, cax=self.ax0)

                cbar.ax.tick_params(labelsize=self.fontsize)

            self.ax0.set_xlim(self.xmin, self.xmax)
            self.ax0.set_ylim(self.ymin, self.ymax)

        # Now plot the cells (possibly on top of the substrate)
        if self.cells_toggle.isChecked():
            if not self.substrates_toggle.isChecked():
                self.fig, (self.ax0) = plt.subplots(1, 1, figsize=(self.figsize_width_svg, self.figsize_height_svg))

            self.svg_frame = frame
            # print('plot_cell_scalar with frame=',self.svg_frame)
            self.plot_cell_scalar(self.svg_frame)

def main():
    output_dir = "output"
    current_frame = 1
    axes_fixed = False
    colorbar_name = "RdBu"
    show_colorbar = False
    scalar_name = "a_i"
    beta_gamma = False
    plot_xmin = 0.0
    plot_xmax = 100.0
    plot_ymin = 0.0
    plot_ymax = 100.0
    try:
        parser = argparse.ArgumentParser(description='Monolayer plot ')

    #   print("Usage:  args = output_dir colorbar_name frame show_colorbar [xmin xmax ymin ymax]")

        parser.add_argument("-o ", "--output_dir", type=str, help="output directory")
        parser.add_argument("-f ", "--frame", type=int, help="current frame to plot")
        # parser.add_argument("-a ", "--axes_fixed", type=bool, help="fixed axes")
        parser.add_argument("-a ", "--axes_fixed", dest="axes_fixed", help="fixed axes", action="store_true")
        parser.add_argument("-c ", "--colorbar_name", type=str, help="mpl colorbar name")
        parser.add_argument("-b ", "--show_colorbar", dest="show_colorbar", help="show colorbar", action="store_true")
        parser.add_argument("-s ", "--scalar_name", type=str, help="scalar value [a_i]")
        parser.add_argument("-x0 ", "--xmin", type=float, help="plot xmin")
        parser.add_argument("-x1 ", "--xmax", type=float, help="plot xmin")
        parser.add_argument("-y0 ", "--ymin", type=float, help="plot ymin")
        parser.add_argument("-y1 ", "--ymax", type=float, help="plot ymax")

        # args = parser.parse_args()
        args, unknown = parser.parse_known_args()
        print("args=",args)
        print("unknown=",unknown)
        if unknown:
            print("len(unknown)= ",len(unknown))
            print("Invalid argument(s): ",unknown)
            print("Use '--help' to see options.")
            sys.exit(-1)
        if args.output_dir:
            output_dir = args.output_dir
        if args.frame:
            current_frame = args.frame
        if args.axes_fixed:
            # axes_fixed = args.axes_fixed
            axes_fixed= True
        if args.colorbar_name:
            colorbar_name = args.colorbar_name
        if args.show_colorbar:
            show_colorbar= True
        if args.scalar_name:
            scalar_name= args.scalar_name
        if args.xmin:
            plot_xmin = args.xmin
            print("args.xmin= ",plot_xmin)
        if args.xmax:
            plot_xmax = args.xmax
            print("args.xmax= ",plot_xmax)
        if args.ymin:
            plot_ymin = args.ymin
            print("args.ymin= ",plot_ymin)
        if args.ymax:
            plot_ymax = args.ymax
            print("args.ymax= ",plot_ymax)

    except:
        print("Error parsing command line args.")
        sys.exit(-1)

    print(f'output_dir={output_dir}, current_frame={current_frame}, axes_fixed={axes_fixed}, colorbar={colorbar_name}, show_colorbar={show_colorbar}, xmax={plot_xmax} ')
    
    # show_colorbar = True
    vis_tab = Vis(output_dir,current_frame,axes_fixed,colorbar_name,show_colorbar,scalar_name, plot_xmin,plot_xmax,plot_ymin,plot_ymax )
	
if __name__ == '__main__':
    main()