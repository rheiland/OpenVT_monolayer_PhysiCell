
import os
#import subprocess

# Roman had beta in X (columns), gamma in Y (rows)
#gamma_vals = [0.5, 0.6, 0.7, 0.8, 0.9]
#gamma_vals = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
gamma_vals = [0.3, 0.5, 0.7, 0.8, 0.9]
# gamma_vals = [0.3, 0.5]
# gamma_vals = [0.3]
print("type(gamma_vals)=", type(gamma_vals))
print(gamma_vals)
gamma_r = list(reversed(gamma_vals))
print("type(gamma_r)=", type(gamma_r))
print(gamma_r)

# os.system("rm results.csv")
os.system("echo 'N,R,A,C,w,g' > results.csv")
os.system("echo 'N,R,A,C,w,g' > results_beta_gamma_params.txt")
# python analysis/get_cells_csv.py run01_nbrs
# analysis/metrics cells.csv >> results.csv

for gamma in gamma_r:
    # for beta in [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]:
    for beta in [0.6, 0.7, 0.8, 0.9, 0.95]:
        output_dir = "out_b" + str(beta) + "_g" + str(gamma)
        cmd = "python analysis/get_cells_csv.py " + output_dir
        os.system(cmd)
        cmd = "analysis/metrics cells.csv >> results.csv"
        os.system(cmd)

        cmd = f"echo {output_dir[4:]} >> results_beta_gamma_params.txt"
        os.system(cmd)

print("\n--> results.csv")
print("--> results_beta_gamma_params.txt")