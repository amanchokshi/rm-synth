"""
Compare a set of RM Spectra
"""

from plot_rm import pogs_obj_loc, read_rm_cube
from matplotlib import pyplot as plt
import matplotlib.patheffects as pe
import matplotlib.pylab as pl
import numpy as np

# Read POGS catalog and create source skycoord object
pogs_pos = pogs_obj_loc("POGSII-EG-321", "../slurm/iono/POGS-II_ExGal.fits")

# List of run_rts --tag names
patches = ["fee", "fee_noflag", "analytic", "analytic_noflag"]

plt.style.use("seaborn")

colors = pl.cm.Spectral(np.linspace(0.05,0.9,len(patches)))

fig = plt.figure(figsize=(9, 6))
ax = fig.add_subplot(1, 1, 1)

for i, patch in enumerate(patches):
    ra_x, dec_y, phi_z, phi, data_p, wcs = read_rm_cube(pogs_pos, "rts_imgr_", f"../data/1086351512/{patch}/run_i/cubes_iono")
    ax.plot(phi, data_p[dec_y, ra_x, :], label=f"{patch}", color=colors[i], alpha=0.8, path_effects=[pe.Stroke(linewidth=2.4, foreground="k"), pe.Normal()])

ax.set_xlabel("Faraday Depth [rad/m$^2$]")
ax.set_ylabel("Polarized Flux Density [Jy/PSF/RMSF]")
ax.set_title("POGSII-EG-321 FEE vs Analytic RM effects")

leg = plt.legend(frameon=True, markerscale=1, handlelength=1)
leg.get_frame().set_facecolor("white")
for le in leg.legendHandles:
    le.set_alpha(1)

plt.tight_layout()
plt.savefig(f"fee_analytic_test.png", dpi=300)
