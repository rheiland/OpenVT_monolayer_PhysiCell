# This script provides simple parameter exploration functionality. The script creates
# a new folder (subdirectory) for each set of parameters, makes changes to a default 
# configuration (.xml) file using specified parameter values (in an accompanying .txt file),
# copies the new config file into the new folder, then
# runs the simulation (in the background) which writes results into the new folder.
# 

import xml.etree.ElementTree as ET
from shutil import copyfile
import os
import sys
import subprocess

# print(len(sys.argv))
if (len(sys.argv) < 2):
  usage_str = "Usage: %s <exec_pgm>" % (sys.argv[0])
  print(usage_str)
  print("e.g.:  python param_sweep.py project")
  exit(1)
else:
   exec_pgm = sys.argv[1]

background_str = " &"  # works on Unix
background_str = " "  # works on Unix
if sys.platform == 'win32':
    background_str = ""


# xml_file_in = 'config/monolayer_gamma_beta.xml'
#xml_file_in = 'mono_gamma_beta.xml'

# xml_file_in = 'monolayer_gamma_beta_custom_volume.xml'
# xml_file_out = 'gamma_beta_sims_custom_volume.xml'
# xml_file_in = 'monolayer_gamma_beta_cell_area.xml'
# xml_file_in = 'monolayer_10K_5x_slower.xml'
xml_file_in = 'monolayer_linear_growth_5x5.xml'

xml_file_out = 'gamma_beta_sims_cell_area.xml'
copyfile(xml_file_in, xml_file_out)
tree = ET.parse(xml_file_out)
xml_root = tree.getroot()
# first_time = True
output_dirs = []

# Roman's latest 2/27/26
# 0, 95, 99, 99.9, 99.95 (same percentiles for beta and gamma)

# Beware trailing zeros!

# NEW - 3/19/26
#a_i: x values for %s: 0.49904, 0.97772, 0.98847, 0.99436, 0.99554
#f_i: x values for %s: 0.00000, 0.71114, 0.82560, 0.90293, 0.91095

# --> 3 decimal digits

# do partial initially
# beta_vals  = [0.0, 0.978, 0.988, 0.994, 0.996]   # use 0.0 instead of 0.499 ??
# gamma_vals = [0.0, 0.711, 0.826, 0.903, 0.911]

# NEW: 3/27/26
# (base) M1P~/git/PhysiCell_monolayer/PhysiCell_mech_grid_xml$ python beta/chosen_CDF_pct.py 100 cdf_1000cells_linear_growth  bg00_cells1000_  f_i
# -->
# x values for %s: 0.00000, 0.72514, 0.83539, 0.90324, 0.91286
gamma_vals = [0.0, 0.725, 0.835, 0.903, 0.913]

# and:
# (base) M1P~/git/PhysiCell_monolayer/PhysiCell_mech_grid_xml$ python beta/chosen_CDF_pct.py 100 cdf_1000cells_linear_growth  bg00_cells1000_  a_i
# -->
# x values for %s: 0.51875, 0.98036, 0.98991, 0.99496, 0.99631
beta_vals  = [0.0, 0.980, 0.990, 0.995, 0.996]   # use 0.0 instead of 0.xxx



# gamma_vals = [0.0, 0.711, 0.826]
# gamma_vals = [ 0.903, 0.911]

# beta_vals  = [0.996]   # use 0.0 instead of 0.499 ??
# gamma_vals = [0.97]
# gamma_vals = [0.95]

for beta in beta_vals:

# for beta in [0.09078]:
# for beta in [0.95]:
    # for gamma in [0.2, 0.5, 0.93, 0.95, 0.96]:
    # for gamma in [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]:

    # for gamma in [0.3, 0.5, 0.7, 0.8, 0.9]:
    # for gamma in [0.0, 0.68266, 0.77982, 0.81597, 0.86907]:
    # for gamma in [0.9081]:
    for gamma in gamma_vals:
        # if beta==0.0 and gamma==0.0:
        #     continue
        # if beta==0.0 and gamma==0.711:
        #     continue
    # for gamma in [0.00000]:
        # folder_name = "out_b" + str(beta) + "_g" + str(gamma)
        # folder_name = "out_new_b" + str(beta) + "_g" + str(gamma)
        folder_name = "out_cell_area_b" + str(beta) + "_g" + str(gamma)
        output_dirs.append(folder_name)
        if (not os.path.exists(folder_name)):
            print("--- mkdir ", folder_name)
            os.makedirs(folder_name)

        xml_file_out = os.path.join(folder_name, 'config.xml')  # copy config file into the output dir

        print('---write config file (and start sim): ', xml_file_out)
        # tree.write(xml_file_out)   # will create folder_name/config.xml
        log_file = folder_name + ".log"  
        cmd =  exec_pgm + " " + xml_file_out + " > " + log_file + " " + background_str
        # print("----- cmd = ",cmd)
        # os.system(cmd)   # <------ Execute the simulation
        # subprocess.Popen([exec_pgm, xml_file_out])
        # with open(log_file,"w") as outf:
        #     subprocess.Popen([exec_pgm, xml_file_out],stdout=outf)

        try:
            xml_root.find('.//' + 'folder').text = str(folder_name)   # beware of rules folder!
        except:
            print("--- Error setting output folder")
            exit(-1)

        try:
            xml_root.find('.//' + 'beta_threshold').text = str(beta)
            xml_root.find('.//' + 'gamma_threshold').text = str(gamma)
        except:
            print("--- Error setting beta or gamma")
            exit(-1)

        xml_file_out = os.path.join(folder_name, 'config.xml')  # copy config file into the output dir
        tree.write(xml_file_out)   # will create folder_name/config.xml

        print("----- cmd = ",cmd)
        os.system(cmd)   # <------ Execute the simulation

        # exit(-1)

if False:
  with open(params_file) as f:
    for line in f:
        # print(len(line),line)
        print(line, end="")
        if (line[0] == '#'):
            continue
        (key, val) = line.split()
        if (key == 'run_it'):
            # write the config file to the previous folder (output) dir and start a simulation
            # print('---write config file and start its sim')
            print('---write config file (and start sim): ', xml_file_out)
            tree.write(xml_file_out)   # will create folder_name/config.xml
            log_file = folder_name + ".log"  
            cmd =  exec_pgm + " " + xml_file_out + " > " + log_file + " " + background_str
            print("----- cmd = ",cmd)
            os.system(cmd)   # <------ Execute the simulation
            # subprocess.Popen([exec_pgm, xml_file_out])
            # with open(log_file,"w") as outf:
            #     subprocess.Popen([exec_pgm, xml_file_out],stdout=outf)
        elif ('.' in key):
            k = key.split('.')
            uep = xml_root
            for idx in range(len(k)):
                uep = uep.find('.//' + k[idx])  # unique entry point (uep) into xml
#                print(k[idx])
            uep.text = val
        else:
            if (key == 'folder'):
                folder_name = val
                output_dirs.append(folder_name)
                if (not os.path.exists(folder_name)):
                    print("--- parsed 'folder', makedir " + folder_name)
                    os.makedirs(folder_name)
                # xml_file_out = folder_name + '/config.xml'  # copy config file into the output dir
                xml_file_out = os.path.join(folder_name, 'config.xml')  # copy config file into the output dir

            try:
                xml_root.find('.//' + key).text = val
            except:
                print("--- Error: could not find ",key," in .xml\n")
                sys.exit(1)

print("\n ------\n Your output results will appear in these directories:\n   ",output_dirs)
print("and check for a .log file of each name for your terminal output from each simulation.\n")
