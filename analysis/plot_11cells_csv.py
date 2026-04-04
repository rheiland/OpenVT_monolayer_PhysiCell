
import numpy as np
import matplotlib.pyplot as plt
import sys

CSV_FILE = sys.argv[1] if len(sys.argv) > 1 else "data.csv"

data = np.genfromtxt(CSV_FILE, delimiter=",", skip_header=1)
time  = data[:, 0]
width = data[:, 1]

fig, ax = plt.subplots()
ax.plot(time, width)
ax.set_ylim(5, 10)
ax.set_xlabel("time")
ax.set_ylabel("width")
ax.axvline(1, linestyle="--", color="blue")
ax.axhline(9, linestyle="--", color="blue")
plt.tight_layout()
# plt.savefig(CSV_FILE.replace(".csv", "_plot.png"), dpi=150)
plt.show()
