"""
An example from Brentjens & Bruyn 2005.

Attempting to recreate Fig. 1 which is a toy model
of a polarized source at 𝛟1 = +10 rad/m^2 with
polarized flux 0.25 Jy/beam. Uniform slab of
galactic foreground at 𝛟fg = 2 rad/m^2
"""

import numpy as np
from matplotlib import pyplot as plt

# Faraday depth of foreground and polarized source
phi_fg = 2
phi_1 = 10

# Lambda squared array
lam2 = np.linspace(0.01, 1, 100)

fg = (1 / (2 * phi_fg * lam2)) * np.sin(2 * phi_fg * lam2)
s_re = 0.25 * np.cos(2 * phi_1 * lam2)
s_im = 0.25 * np.sin(2 * phi_1 * lam2)
source = s_re + s_im * 1j

los = source + fg

uq = np.rad2deg(np.arctan(los.imag / los.real))

uq_min = np.amin(uq)

#  chi = 1 * (((uq - uq_min) % 90) - 0)
chi = 0.5 * uq

print(chi)


# Plotting stuff
plt.style.use("seaborn")

fig = plt.figure(figsize=(9, 6))
ax = fig.add_subplot(1, 1, 1)

#  ax.plot(lam2, fg, color="seagreen", label="P$_{fg}$")
#  ax.plot(lam2, np.abs(los), color="indigo", label="||P$_{1}$||")
ax.plot(lam2, chi, color="indigo", label="$\chi$")

ax.set_xlabel("$\lambda^{2} \\ [m^{2}]$")
#  ax.set_ylabel("Flux [Jy]")
ax.set_ylabel("$\chi$ [deg]")
leg = plt.legend(frameon=True, markerscale=1, handlelength=1)
leg.get_frame().set_facecolor("white")
for le in leg.legendHandles:
    le.set_alpha(1)

plt.tight_layout()
plt.show()
