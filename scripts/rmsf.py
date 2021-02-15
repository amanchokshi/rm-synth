import numpy as np
from matplotlib import pyplot as plt

file = "./rm_cubes/wsclean_dirty_rmsf.txt"
data = np.loadtxt(file)

phi = data[:, 0]
q = data[:, 1]
u = data[:, 2]
p = data[:, 3]


plt.style.use("seaborn")
plt.plot(phi, p, label=r"$ \vert R \vert $")
plt.plot(phi, q, label=r"$ real(R) $")
plt.plot(phi, u, label=r"$ imag(R) $")
#  plt.xlim([-50, 50])
plt.xlabel("Faraday Depth")
plt.ylabel("RMSF")
plt.legend()
plt.tight_layout()
plt.show()
