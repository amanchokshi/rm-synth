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
    #  ax.set_ylim([0.0, 7])
    ax.set_xlabel("Dipole Gain Amplitude [0, 1]")
    #  ax.set_yticks([])

    leg = ax.legend(loc="upper left", frameon=True, markerscale=4, handlelength=1)
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(0.9)


def log_chisq(model, data, mask):
    model[mask] = np.nan
    data[mask] = np.nan

    return np.log(np.nansum(np.square(model - data)))


if __name__ == "__main__":

    from pathlib import Path

    #  beam_min = np.load("../data/mcmc/S06XX_beam_min_1024_walk.npy")
    tiles = [
        "S06XX",
        "S06YY",
        "S07XX",
        "S07YY",
        "S08XX",
        "S08YY",
        "S09XX",
        "S09YY",
        "S10XX",
        "S10YY",
        "S12XX",
        "S12YY",
        "S29XX",
        "S29YY",
        "S30XX",
        "S30YY",
        "S31XX",
        "S31YY",
        "S32XX",
        "S32YY",
        "S33XX",
        "S33YY",
        "S34XX",
        "S34YY",
        "S35XX",
        "S35YY",
        "S36XX",
        "S36YY",
    ]

    chi_sq_tiles = []

    for tile in tiles:

        out_dir = Path(f"../data/mcmc/beam_min_masked_512/{tile}/")
        out_dir.mkdir(parents=True, exist_ok=True)

        beam_min = np.load(
            f"../data/mcmc/beam_min_masked_512/data/{tile}_beam_min_512_walk_mask.npy"
        )

        dipoles = [
            "$A_0$",
            "$A_1$",
            "$A_2$",
            "$A_3$",
            "$A_4$",
            "$A_5$",
            "$A_6$",
            "$A_7$",
            "$A_8$",
            "$A_9$",
            "$A_{10}$",
            "$A_{11}$",
            "$A_{12}$",
            "$A_{13}$",
            "$A_{14}$",
            "$A_{15}$",
        ]

        plt.style.use("seaborn")

        plt.rcParams.update(
            {
                "font.size": 16,
                "font.family": "serif",
                "font.serif": "Times New Roman",
                "text.usetex": True,
            }
        )
        fig, axs = plt.subplots(4, 4, sharex=True, sharey=True, figsize=(16, 9))
        plot_beam_amp(beam_min, dipoles, 0, fig, axs[0, 0])
        plot_beam_amp(beam_min, dipoles, 1, fig, axs[0, 1])
        plot_beam_amp(beam_min, dipoles, 2, fig, axs[0, 2])
        plot_beam_amp(beam_min, dipoles, 3, fig, axs[0, 3])
        plot_beam_amp(beam_min, dipoles, 4, fig, axs[1, 0])
        plot_beam_amp(beam_min, dipoles, 5, fig, axs[1, 1])
        plot_beam_amp(beam_min, dipoles, 6, fig, axs[1, 2])
        plot_beam_amp(beam_min, dipoles, 7, fig, axs[1, 3])
        plot_beam_amp(beam_min, dipoles, 8, fig, axs[2, 0])
        plot_beam_amp(beam_min, dipoles, 9, fig, axs[2, 1])
        plot_beam_amp(beam_min, dipoles, 10, fig, axs[2, 2])
        plot_beam_amp(beam_min, dipoles, 11, fig, axs[2, 3])
        plot_beam_amp(beam_min, dipoles, 12, fig, axs[3, 0])
        plot_beam_amp(beam_min, dipoles, 13, fig, axs[3, 1])
        plot_beam_amp(beam_min, dipoles, 14, fig, axs[3, 2])
        plot_beam_amp(beam_min, dipoles, 15, fig, axs[3, 3])

        # Only outer labels
        for ax in fig.get_axes():
            ax.label_outer()

        plt.tight_layout()
        plt.savefig(f"{out_dir}/{tile}_beam_min.png")
        plt.close()

        amps_sat = np.median(beam_min, axis=0)

        # Hyperbeam settings
        nside = 32
        freq = 138e6
        delays = [0] * 16
        norm_to_zenith = True

        # Zenith angle and Azimuth of healpix pixels
        za, az = bu.healpix_za_az(nside=nside)

        # Make a new beam object
        beam = mwa_hyperbeam.FEEBeam()

        if "XX" in tile:
            pol = "XX"
        else:
            pol = "YY"

        # Load satellite beam map
        data_sat = np.load(f"../data/embers_healpix/{tile}_rf1{pol}_0.npz")["beam_map"][
            : az.shape[0]
        ]

        # Create mask based on -30dB threshold of perfect FEE model
        jones_perfect = beam.calc_jones_array(
            az, za, freq, delays, [1.0] * 16, norm_to_zenith
        )
        unpol_perfect = bu.makeUnpolInstrumentalResponse(jones_perfect, jones_perfect)

        if pol == "XX":
            fee_perfect = 10 * np.log10(np.real(unpol_perfect[:, 0]))
            mask_30dB = np.where(fee_perfect < -30)
        else:
            fee_perfect = 10 * np.log10(np.real(unpol_perfect[:, 3]))
            mask_30dB = np.where(fee_perfect < -30)

        jones_min = beam.calc_jones_array(
            az, za, freq, delays, amps_sat, norm_to_zenith
        )
        unpol_min = bu.makeUnpolInstrumentalResponse(jones_min, jones_min)

        if pol == "XX":
            fee_min = 10 * np.log10(np.real(unpol_min[:, 0]))
        else:
            fee_min = 10 * np.log10(np.real(unpol_min[:, 3]))

        fee_chisq = log_chisq(fee_perfect, data_sat, mask_30dB)
        fee_min_chisq = log_chisq(fee_min, data_sat, mask_30dB)

        chi_sq_tiles.append([fee_chisq, fee_min_chisq])

    # Hyperbeam settings
    nside = 32
    freq = 138e6
    delays = [0] * 16
    norm_to_zenith = True

    # Zenith angle and Azimuth of healpix pixels
    za, az = bu.healpix_za_az(nside=nside)

    # Make a new beam object
    beam = mwa_hyperbeam.FEEBeam()

    # Create mask based on -30dB threshold of perfect FEE model
    jones_perfect = beam.calc_jones_array(
        az, za, freq, delays, [1.0] * 16, norm_to_zenith
    )
    unpol_perfect = bu.makeUnpolInstrumentalResponse(jones_perfect, jones_perfect)

    model_perfect_XX = 10 * np.log10(np.real(unpol_perfect[:, 0]))
    mask_30dB_XX = np.where(model_perfect_XX < -30)

    amps_flagged = [0.0] + [1.0] * 15
    log_chi_sq = []
    for i in range(16):
        amps_f = np.roll(amps_flagged, i)

        jones_min = beam.calc_jones_array(az, za, freq, delays, amps_f, norm_to_zenith)
        unpol_min = bu.makeUnpolInstrumentalResponse(jones_min, jones_min)

        data_flagged_XX = 10 * np.log10(np.real(unpol_min[:, 0]))

        lchs = log_chisq(model_perfect_XX, data_flagged_XX, mask_30dB_XX)

        log_chi_sq.append(lchs)

    chi_sq_tiles = np.array(chi_sq_tiles)
    plt.style.use("seaborn")
    plt.rcParams.update(
        {
            "font.size": 16,
            "font.family": "serif",
            "font.serif": "Times New Roman",
            "text.usetex": True,
        }
    )
    fig, ax = plt.subplots(1, 1, figsize=(11, 7))
    ax.plot(range(28), chi_sq_tiles[:, 0], "-D", color="#008891", label="FEE Perfect")
    ax.plot(range(28), chi_sq_tiles[:, 1], "-s", color="#00587a", label="FEE Minimized")
    ax.hlines(
        np.mean(log_chi_sq),
        0,
        27,
        colors="#687980",
        linestyles="dashed",
        label="Avg Flagged Dipole",
    )
    ax.set_xticks(np.arange(28))
    ax.set_xticklabels(tiles)
    ax.tick_params(axis="x", rotation=90)
    ax.set_ylabel("Log Chi Square")
    ax.set_title("MWA Satellite Beam Map Gain Minimization : Log Chi Square Test")

    leg = ax.legend(loc="upper right", frameon=True)
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    plt.tight_layout()
    plt.savefig("../data/mcmc/beam_min_masked_512/MWA_Min_Log_Chisq.png")
    #  plt.show()

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
    #         "font.family": "serif",
    #         "font.serif": "Times New Roman",
    #         "text.usetex": True,
    #     }
    # )
    # fig = plt.subplots(1, 1, figsize=(11, 7))
    # plt.plot(range(16), log_chi_sq, "-o", color="midnightblue")
    # plt.tight_layout()
    # plt.savefig("../data/mcmc/beam_min_masked_512/Flagged_FEE_Log_Chisq.png")
