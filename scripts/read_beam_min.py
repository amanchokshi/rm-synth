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
                #  "xtick.major.size": 12,
                #  "ytick.major.size": 12,
                #  "xtick.major.width": 1,
                #  "ytick.major.width": 1,
                #  "ytick.minor.size": 5,
                #  "xtick.minor.size": 5,
                #  "axes.linewidth": 1,
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
        data_sat = np.load(f"../data/embers_healpix/{tile}_rf1{pol}_0.npz")["beam_map"]

        # Create mask based on -30dB threshold of perfect FEE model
        jones_perfect = beam.calc_jones_array(
            az, za, freq, delays, [1.0] * 16, norm_to_zenith
        )
        unpol_perfect = bu.makeUnpolInstrumentalResponse(jones_perfect, jones_perfect)

        if pol == "XX":
            model_perfect = 10 * np.log10(np.real(unpol_perfect[:, 0]))
            mask_30dB = np.where(model_perfect < -30)
        else:
            model_perfect = 10 * np.log10(np.real(unpol_perfect[:, 3]))
            mask_30dB = np.where(model_perfect < -30)

        fee_perfect = np.zeros(hp.nside2npix(nside))
        fee_perfect[: model_perfect.shape[0]] = model_perfect

        jones_min = beam.calc_jones_array(
            az, za, freq, delays, amps_sat, norm_to_zenith
        )
        unpol_min = bu.makeUnpolInstrumentalResponse(jones_min, jones_min)

        if pol == "XX":
            data_min = 10 * np.log10(np.real(unpol_min[:, 0]))
        else:
            data_min = 10 * np.log10(np.real(unpol_min[:, 3]))

        fee_min = np.zeros(hp.nside2npix(nside))
        fee_min[: data_min.shape[0]] = data_min

        plt.rcParams.update(plt.rcParamsDefault)
        fig = plt.subplots(1, 1, figsize=(7, 8))
        bu.plot_healpix(
            data_map=fee_perfect, vmin=-50, vmax=0, cmap="viridis", title="FEE Perfect",
        )
        plt.tight_layout()
        plt.savefig(f"{out_dir}/FEE.png")
        plt.close()

        fig = plt.subplots(1, 1, figsize=(7, 8))
        bu.plot_healpix(
            data_map=fee_min, vmin=-50, vmax=0, cmap="viridis", title="FEE Minimized",
        )
        plt.tight_layout()
        plt.savefig(f"{out_dir}/FEE_Min.png")
        plt.close()

        fig = plt.subplots(1, 1, figsize=(7, 8))
        bu.plot_healpix(
            data_map=data_sat,
            vmin=-50,
            vmax=0,
            cmap="viridis",
            title=f"{tile} Satellite Map",
        )
        plt.tight_layout()
        plt.savefig(f"{out_dir}/{tile}_Sat.png")
        plt.close()

        fig = plt.subplots(1, 1, figsize=(7, 8))
        data_res = fee_perfect - data_sat
        fp = fee_perfect
        ds = data_sat
        fp[mask_30dB] = np.nan
        ds[mask_30dB] = np.nan
        log_chi_sq = np.log(np.nansum(np.square(fp - ds)))
        data_res[mask_30dB] = np.nan
        bu.plot_healpix(
            data_map=data_res,
            vmin=-5,
            vmax=5,
            cmap="RdYlGn",
            title=f"FEE - {tile} Satellite Map : Log Chisq = {log_chi_sq:.3f}",
        )
        plt.tight_layout()
        plt.savefig(f"{out_dir}/FEE-{tile}_Sat.png")
        plt.close()

        fig = plt.subplots(1, 1, figsize=(7, 8))
        min_res = fee_min - data_sat
        min_res[mask_30dB] = np.nan
        fm = fee_min
        ds = data_sat
        fm[mask_30dB] = np.nan
        ds[mask_30dB] = np.nan
        log_chi_sq = np.log(np.nansum(np.square(fm - ds)))
        bu.plot_healpix(
            data_map=min_res,
            vmin=-5,
            vmax=5,
            cmap="RdYlGn",
            title=f"FEE Min - {tile} Satellite Map : Log Chisq = {log_chi_sq:.3f}",
        )
        plt.tight_layout()
        plt.savefig(f"{out_dir}/FEE_Min-{tile}_Sat.png")
        plt.close()
