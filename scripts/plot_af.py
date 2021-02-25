from matplotlib import pyplot as plt

from plot_rm import pogs_obj_loc, read_noise, read_rm_cube

pogs_pos = pogs_obj_loc("POGSII-EG-321", "../slurm/iono/POGS-II_ExGal.fits")

ra_a, dec_a, phi_a, phi_a, data_a, wcs_a = read_rm_cube(
    pogs_pos, "ana_", "../data/noise_321"
)
n_a = read_noise(pogs_pos, "ana_", "../data/noise_321")
ra_f, dec_f, phi_f, phi_f, data_f, wcs_f = read_rm_cube(
    pogs_pos, "fee_", "../data/noise_321"
)
n_f = read_noise(pogs_pos, "fee_", "../data/noise_321")

plt.style.use("seaborn")

fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(1, 1, 1)
ax.plot(
    phi_f,
    data_f[dec_f, ra_f, :] - n_f,
    label="FEE",
    linewidth=2,
    color="#28527a",
    alpha=0.7,
)
ax.plot(
    phi_a,
    data_a[dec_a, ra_a, :] - n_a,
    label="ANA",
    linewidth=1,
    linestyle="dashed",
    dash_capstyle="round",
    color="#ef4f4f",
)
ax.set_xlabel("Faraday Depth [rad/m$^2$]")
ax.set_ylabel("Polarized Flux Density [Jy/PSF/RMSF]")
ax.set_title("POGSII-EG-321 Analytic vs FEE Beam")
ax.grid(True, color="white", linewidth=1.2, alpha=0.9, ls="dotted")

leg = plt.legend(frameon=True, markerscale=1, handlelength=1)
leg.get_frame().set_edgecolor('#22222233')
leg.get_frame().set_facecolor("none")
for le in leg.legendHandles:
    le.set_alpha(1)

plt.tight_layout()
plt.savefig("rm_spec_ana_fee.png", transparent=True, dpi=300)
