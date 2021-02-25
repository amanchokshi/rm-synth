"""
Compare a set of RM Spectra
"""

import numpy as np
from matplotlib import patheffects as pe
from matplotlib import pylab as pl
from matplotlib import pyplot as plt

from plot_rm import pogs_obj_loc, read_noise, read_rm_cube

# Read POGS catalog and create source skycoord object
pogs_pos = pogs_obj_loc("POGSII-EG-321", "../slurm/iono/POGS-II_ExGal.fits")
obsid = "1120300232"

# List of run_rts --tag names
#  tags = ["fee_p321", "fee_p321_nodflag", "ana_p321", "ana_p321_nodflag"]
#  tags = ["fee_p321", "ana_p321"]
tags = ["ana_p321_nodflag", "ana_p321"]
#  tags = ["fee_p321", "fee_p321_nodflag"]

plt.style.use("seaborn")

#  colors = pl.cm.Spectral(np.linspace(0.05, 0.9, len(tags)))
#  colors = pl.cm.Spectral(np.linspace(0.02, 0.78, 2))
#  colors = ["#cc561e", "#008891"] # ANA vs FEE
colors = ["#184d47", "#96bb7c"]  # ANA vs Dipole
#  colors = ["#583d72", "#de4463"]  # FEE vs Dipole

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1)

for i, tag in enumerate(tags):
    ra_x, dec_y, phi_z, phi, data_p, wcs = read_rm_cube(
        pogs_pos,
        "rts_imgr_",
        f"../data/{obsid}/{tag}"
        #  pogs_pos, "rts_imgr_", f"../data/{obsid}/{tag}/imgs/cubes_iono"
    )
    noise = read_noise(pogs_pos, "", f"../data/{obsid}/{tag}")
    ax.plot(
        phi,
        data_p[dec_y, ra_x, :] - noise,
        label=f"{tag}",
        #  color=colors[::-2][i],
        linewidth=2,
        color=colors[i],
        alpha=0.9,
        #  path_effects=[pe.Stroke(linewidth=2.4, foreground="k"), pe.Normal()],
    )

ax.set_xlabel("Faraday Depth [rad/m$^2$]")
ax.set_ylabel("Polarized Flux Density [Jy/PSF/RMSF]")
ax.set_title("POGSII-EG-321 Analytic vs Dipole Flag RM effects")
#  ax.set_title("POGSII-EG-321 Analytic vs FEE RM Spectra")

ax.set_ylim([-0.013, 0.13])
#  ax.set(ylabel=None)
#  ax.set(yticklabels=[])

ax.grid(True, color="white", linewidth=1.2, alpha=0.9, ls="dotted")

leg = plt.legend(frameon=True, markerscale=1, handlelength=1)
leg.get_frame().set_edgecolor("#22222233")
leg.get_frame().set_facecolor("none")
for le in leg.legendHandles:
    le.set_alpha(1)

plt.tight_layout()
plt.savefig("ana_dipole.png", dpi=300, transparent=True)
