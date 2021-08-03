import json

import healpy as hp
import mwa_hyperbeam
import numpy as np
from matplotlib import pyplot as plt

import beam_utils as bu
from beam_min_combinations import write_json
from beam_minimize import beam_mask


def log_chisq(model, data, mask):
    model[mask] = np.nan
    data[mask] = np.nan

    return np.log(np.nansum(np.square(model - data)))


if __name__ == "__main__":

    from pathlib import Path

    tiles = [
        "S06XX_rf1XX",
        "S06YY_rf1YY",
        "S07XX_rf1XX",
        "S07YY_rf1YY",
        "S08XX_rf1XX",
        "S08YY_rf1YY",
        "S09XX_rf1XX",
        "S09YY_rf1YY",
        "S10XX_rf1XX",
        "S10YY_rf1YY",
        "S12XX_rf1XX",
        "S12YY_rf1YY",
        "S29XX_rf1XX",
        "S29YY_rf1YY",
        "S30XX_rf1XX",
        "S30YY_rf1YY",
        "S31XX_rf1XX",
        "S31YY_rf1YY",
        "S32XX_rf1XX",
        "S32YY_rf1YY",
        "S33XX_rf1XX",
        "S33YY_rf1YY",
        "S34XX_rf1XX",
        "S34YY_rf1YY",
        "S35XX_rf1XX",
        "S35YY_rf1YY",
        "S36XX_rf1XX",
        "S36YY_rf1YY",
    ]

    # Each element of this list contains the set of chisq values of
    # the satellite map compared to the full FEE beam, and the minimized (median)
    # FEE beam
    chi_sq_tiles = []

    # Chisqs from all possible dipole peak combinations
    chi_sq_comb_min = []
    amp_comb_min = {}

    for tile in tiles:

        out_dir = Path(f"../data/beam_min_1024_masked/{tile}/")
        out_dir.mkdir(parents=True, exist_ok=True)

        # All 1024 min dipole gains from beam_minimize.py, started from random initial conditions
        # Histogram this data to get an idea of prefered dipole gains
        beam_min = np.load(
            f"../data/beam_min_1024_masked/raw/{tile}_beam_min_1024_walk_mask.npy"
        )

        # KDE fit to above histogram, and upto two highest peaks selected per dipole
        # Chisq of ll combinations of possible peaks evaluated and saved to json dict
        with open(
            f"../data/beam_min_1024_masked/raw/{tile}_amp_combinations.json", "r"
        ) as file:
            data = file.read()

            amps_comb = np.array(list(json.loads(data).values()), dtype="object")
            chisq = amps_comb[:, -1]
            chisq_min = min(chisq)

            min_ind = np.where(chisq == chisq_min)[0][0]

            chi_sq_comb_min.append(chisq_min)
            amp_comb_min[f"{tile.split('_')[0]}"] = amps_comb[:, 0][min_ind]

        # Median dipole gains from all 1024 random walkers for beam_minimize.py
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
        map_med = np.load(f"../data/embers_maps/rf1_med_maps/{tile}_med.npy")
        map_mad = np.load(f"../data/embers_maps/rf1_mad_maps/{tile}_mad.npy")

        mask = beam_mask(map_med, map_mad, pol=pol)

        jones_perfect = beam.calc_jones_array(
            az, za, freq, delays, [1.0] * 16, norm_to_zenith
        )
        unpol_perfect = bu.makeUnpolInstrumentalResponse(jones_perfect, jones_perfect)

        if pol == "XX":
            fee_perfect = 10 * np.log10(np.real(unpol_perfect[:, 0]))
        else:
            fee_perfect = 10 * np.log10(np.real(unpol_perfect[:, 3]))

        fee = np.zeros(hp.nside2npix(nside))
        fee[: fee_perfect.shape[0]] = fee_perfect

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

        fee[mask] = np.nan
        fee_min[mask] = np.nan
        map_med[mask] = np.nan
        map_mad[mask] = np.nan

        log_chisq_fee = np.log(np.nansum(np.square(fee - map_med) / map_mad))
        log_chisq_fee_min = np.log(np.nansum(np.square(fee_min - map_med) / map_mad))
        chi_sq_tiles.append([log_chisq_fee, log_chisq_fee_min])


    # # Hyperbeam settings
    # nside = 32
    # freq = 138e6
    # delays = [0] * 16
    # norm_to_zenith = True

    # # Zenith angle and Azimuth of healpix pixels
    # za, az = bu.healpix_za_az(nside=nside)

    # # Make a new beam object
    # beam = mwa_hyperbeam.FEEBeam()

    # # Create mask based on -30dB threshold of perfect FEE model
    # jones_perfect = beam.calc_jones_array(
    #     az, za, freq, delays, [1.0] * 16, norm_to_zenith
    # )
    # unpol_perfect = bu.makeUnpolInstrumentalResponse(jones_perfect, jones_perfect)

    # model_perfect_XX = 10 * np.log10(np.real(unpol_perfect[:, 0]))
    # mask_30dB_XX = np.where(model_perfect_XX < -30)

    # amps_flagged = [0.0] + [1.0] * 15
    # log_chi_sq = []
    # for i in range(16):
    #     amps_f = np.roll(amps_flagged, i)

    #     jones_min = beam.calc_jones_array(az, za, freq, delays, amps_f, norm_to_zenith)
    #     unpol_min = bu.makeUnpolInstrumentalResponse(jones_min, jones_min)

    #     data_flagged_XX = 10 * np.log10(np.real(unpol_min[:, 0]))

    #     lchs = log_chisq(model_perfect_XX, data_flagged_XX, mask_30dB_XX)

    #     log_chi_sq.append(lchs)

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
    ax.plot(
        range(len(tiles)),
        chi_sq_tiles[:, 0],
        "-D",
        color="#008891",
        label="FEE Perfect",
    )
    # ax.plot(
    #     range(len(tiles)),
    #     chi_sq_tiles[:, 1],
    #     "-s",
    #     color="#AA2B1D",
    #     label="FEE Min Median",
    # )
    ax.plot(
        range(len(chi_sq_comb_min)),
        chi_sq_comb_min,
        "-s",
        color="#00587a",
        label="FEE Comb Min",
    )

    # ax.hlines(
    #     np.mean(log_chi_sq),
    #     0,
    #     len(tiles) - 1,
    #     colors="#687980",
    #     linestyles="dashed",
    #     label="Avg Flagged Dipole",
    # )

    tile_names = [
        "$S06XX$",
        "$S06YY$",
        "$S07XX$",
        "$S07YY$",
        "$S08XX$",
        "$S08YY$",
        "$S09XX$",
        "$S09YY$",
        "$S10XX$",
        "$S10YY$",
        "$S12XX$",
        "$S12YY$",
        "$S29XX$",
        "$S29YY$",
        "$S30XX$",
        "$S30YY$",
        "$S31XX$",
        "$S31YY$",
        "$S32XX$",
        "$S32YY$",
        "$S33XX$",
        "$S33YY$",
        "$S34XX$",
        "$S34YY$",
        "$S35XX$",
        "$S35YY$",
        "$S36XX$",
        "$S36YY$",
    ]

    ax.set_xticks(range(len(tile_names)))
    ax.set_xticklabels(tile_names)
    ax.tick_params(axis="x", rotation=90)
    ax.set_ylabel("Log Chi Square")
    ax.set_title("MWA Satellite Beam Map Gain Minimization : Log Chi Square Test")

    leg = ax.legend(loc="upper right", frameon=True)
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    plt.tight_layout()
    #  plt.show()
    plt.savefig("../data/beam_min_1024_masked/MWA_Min_Log_Chisq.png")

    write_json(amp_comb_min, filename="sat_dipole_amps.json", out_dir="../data/beam_min_1024_masked")

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
    # plt.plot(range(len(dipoles)), log_chi_sq, "-o", color="midnightblue")
    # plt.tight_layout()
    # plt.savefig("../data/beam_min_1024_masked/Flagged_FEE_Log_Chisq.png")
