import numpy as np
from matplotlib import pyplot as plt

file = "../data/rts_imgr_rmsf.txt"
data = np.loadtxt(file)

phi = data[:, 0]
q = data[:, 1]
u = data[:, 2]
p = data[:, 3]


plt.style.use("seaborn")
plt.plot(phi, p, color="#207561", linewidth=1.2, label=r"$ \vert R \vert $", zorder=3)
plt.plot(
    phi,
    q,
    alpha=0.7,
    color="#da4302",
    linewidth=2.8207561,
    linestyle=(0, (0.1, 2)),
    dash_capstyle="round",
    label=r"$ real(R) $",
    zorder=2
)
plt.plot(
    phi,
    u,
    alpha=0.93,
    color="#552244",
    linewidth=1.6,
    linestyle="dashed",
    dash_capstyle="round",
    label=r"$ imag(R) $",
    zorder=1
)
plt.xlabel("Faraday Depth")
plt.ylabel("RMSF")
plt.legend()
plt.tight_layout()
#  plt.savefig("rmsf.png")
plt.show()
