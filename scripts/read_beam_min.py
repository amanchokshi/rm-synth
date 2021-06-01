import healpy as hp
import mwa_hyperbeam
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

import beam_utils as bu


def plot_beam_amp(beam_min, dipoles, index, fig, ax):

    cols = plt.cm.Spectral(np.linspace(0, 1, 16))
    kde = stats.gaussian_kde(beam_min[:, index])
    ax.hist(
        beam_min[:, index],
        bins=np.arange(0.4, 1.02, 0.02),
        ec="black",
        density=True,
        color=cols[index],
        alpha=0.99,
        label=f"{dipoles[index]} : $\mu$={np.median(beam_min[:, index]):.2f}",
    )
    ax.plot(
        np.linspace(min(beam_min[:, index]), max(beam_min[:, index]), 128),
        kde(np.linspace(min(beam_min[:, index]), max(beam_min[:, index]), 128)),
        linewidth=2,
    )
    ax.vlines(
        np.median(beam_min[:, index]),
        0,
        kde(np.median(beam_min[:, index])),
        color="black",
        linewidth=2,
        linestyle="dashed",
    )
    ax.set_xlim([0.4, 1.02])
    ax.set_xlabel("Dipole Gain Amplitude [0, 1]")
    ax.set_yticks([])

    leg = ax.legend(loc="upper left", frameon=True, markerscale=4, handlelength=1)
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(0.9)


if __name__ == "__main__":

    beam_min = np.load("../data/mcmc/S06XX_beam_min_1024_walk.npy")

    # dipoles = [
    #     "$A_0$",
    #     "$A_1$",
    #     "$A_2$",
    #     "$A_3$",
    #     "$A_4$",
    #     "$A_5$",
    #     "$A_6$",
    #     "$A_7$",
    #     "$A_8$",
    #     "$A_9$",
    #     "$A_{10}$",
    #     "$A_{11}$",
    #     "$A_{12}$",
    #     "$A_{13}$",
    #     "$A_{14}$",
    #     "$A_{15}$",
    # ]

    # plt.style.use("seaborn")

    # plt.rcParams.update(
    #     {
    #         "font.size": 16,
    #         #  "xtick.major.size": 12,
    #         #  "ytick.major.size": 12,
    #         #  "xtick.major.width": 1,
    #         #  "ytick.major.width": 1,
    #         #  "ytick.minor.size": 5,
    #         #  "xtick.minor.size": 5,
    #         #  "axes.linewidth": 1,
    #         "font.family": "serif",
    #         "font.serif": "Times New Roman",
    #         "text.usetex": True,
    #     }
    # )
    # fig, axs = plt.subplots(4, 4, sharex=True, figsize=(16, 9))
    # plot_beam_amp(beam_min, dipoles, 0, fig, axs[0, 0])
    # plot_beam_amp(beam_min, dipoles, 1, fig, axs[0, 1])
    # plot_beam_amp(beam_min, dipoles, 2, fig, axs[0, 2])
    # plot_beam_amp(beam_min, dipoles, 3, fig, axs[0, 3])
    # plot_beam_amp(beam_min, dipoles, 4, fig, axs[1, 0])
    # plot_beam_amp(beam_min, dipoles, 5, fig, axs[1, 1])
    # plot_beam_amp(beam_min, dipoles, 6, fig, axs[1, 2])
    # plot_beam_amp(beam_min, dipoles, 7, fig, axs[1, 3])
    # plot_beam_amp(beam_min, dipoles, 8, fig, axs[2, 0])
    # plot_beam_amp(beam_min, dipoles, 9, fig, axs[2, 1])
    # plot_beam_amp(beam_min, dipoles, 10, fig, axs[2, 2])
    # plot_beam_amp(beam_min, dipoles, 11, fig, axs[2, 3])
    # plot_beam_amp(beam_min, dipoles, 12, fig, axs[3, 0])
    # plot_beam_amp(beam_min, dipoles, 13, fig, axs[3, 1])
    # plot_beam_amp(beam_min, dipoles, 14, fig, axs[3, 2])
    # plot_beam_amp(beam_min, dipoles, 15, fig, axs[3, 3])

    # # Only outer labels
    # for ax in fig.get_axes():
    #     ax.label_outer()

    # plt.tight_layout()
    # plt.savefig("./beam_min/beam_min.png")

    amps_S06XX = np.median(beam_min, axis=0)

    data_S06XX = np.load("../data/embers_healpix/S06XX_rf1XX_0.npz")["beam_map"]
    #  bu.plot_healpix(data_map=data_S06XX, vmin=-50, vmax=0)
    #  plt.show()
    #  plt.close()

    # Healpix map with given nside
    nside = 32

    # Zenith angle and Azimuth of healpix pixels
    za, az = bu.healpix_za_az(nside=nside)

    # Make a new beam object
    beam = mwa_hyperbeam.FEEBeam()

    # Satellite beam map frequency
    freq = 138e6

    # Zenith Pointing
    delays = [0] * 16

    # Normalize maps to zenith
    norm_to_zenith = True

    amps_perfect = [1.0] * 16
    jones_perfect = beam.calc_jones_array(
        az, za, freq, delays, amps_perfect, norm_to_zenith
    )
    unpol_perfect = bu.makeUnpolInstrumentalResponse(jones_perfect, jones_perfect)
    data_perfect_XX = 10 * np.log10(np.real(unpol_perfect[:, 0]))

    data_perfect = np.zeros(hp.nside2npix(nside))
    data_perfect[: data_perfect_XX.shape[0]] = data_perfect_XX

    jones_min = beam.calc_jones_array(az, za, freq, delays, amps_S06XX, norm_to_zenith)
    unpol_min = bu.makeUnpolInstrumentalResponse(jones_min, jones_min)
    data_min_XX = 10 * np.log10(np.real(unpol_min[:, 0]))

    data_min = np.zeros(hp.nside2npix(nside))
    data_min[: data_min_XX.shape[0]] = data_min_XX

    # bu.plot_healpix(
    #     #  data_map=data_perfect - data_S06XX, vmin=-5, vmax=5, cmap="RdYlGn", title="FEE - Sat S06XX"
    #     data_map=data_min - data_S06XX,
    #     vmin=-5,
    #     vmax=5,
    #     cmap="RdYlGn",
    #     #  title=f"FEE- Sat S06XX : Chisq [{np.log(np.nansum(np.square(data_perfect - data_S06XX))):.3f}]",
    #     title=f"FEE Minimized - Sat S06XX",
    # )
    # plt.tight_layout()
    # plt.savefig("./beam_min/FEE_Min-Sat_S06XX.png")
    # plt.close()

    data_perfect[np.where(data_perfect < -30)] = np.nan

    bu.plot_healpix(
        data_map=data_perfect,
        vmin=-50,
        vmax=0,
    )
    #  plt.tight_layout()
    #  plt.savefig("FEE_Min-SatS06XX.png")
    #  plt.close()
    plt.show()
