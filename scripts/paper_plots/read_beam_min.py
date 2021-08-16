import json
import sys

import healpy as hp
import mwa_hyperbeam
import numpy as np
from matplotlib import pyplot as plt

sys.path.append("..")
import beam_utils as bu
from beam_minimize import beam_mask

if __name__ == "__main__":

    from pathlib import Path

    # The best dipole amp combinations
    sat_amps = "../../data/beam_min_1024_masked/sat_dipole_amps.json"

    with open(sat_amps, "r") as file:
        data = file.read()
        dip_amps_min = json.loads(data)

    tiles = [
        # "S06XX_rf1XX",
        "S06YY_rf1YY",
        # "S07XX_rf1XX",
        # "S07YY_rf1YY",
        # "S08XX_rf1XX",
        # "S08YY_rf1YY",
        # "S09XX_rf1XX",
        # "S09YY_rf1YY",
        # "S10XX_rf1XX",
        # "S10YY_rf1YY",
        # "S12XX_rf1XX",
        # "S12YY_rf1YY",
        # "S29XX_rf1XX",
        # "S29YY_rf1YY",
        # "S30XX_rf1XX",
        # "S30YY_rf1YY",
        # "S31XX_rf1XX",
        # "S31YY_rf1YY",
        # "S32XX_rf1XX",
        # "S32YY_rf1YY",
        # "S33XX_rf1XX",
        # "S33YY_rf1YY",
        # "S34XX_rf1XX",
        # "S34YY_rf1YY",
        # "S35XX_rf1XX",
        # "S35YY_rf1YY",
        # "S36XX_rf1XX",
        # "S36YY_rf1YY",
    ]

    for tile in tiles:

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

        amps_sat = dip_amps_min[f"{tile.split('_')[0]}"]

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
        map_med = np.load(f"../../data/embers_maps/rf1_med_maps/{tile}_med.npy")
        map_mad = np.load(f"../../data/embers_maps/rf1_mad_maps/{tile}_mad.npy")

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

        plt.rcParams.update(plt.rcParamsDefault)
        fig = plt.subplots(1, 1, figsize=(7, 8))
        bu.plot_healpix(
            data_map=fee, vmin=-50, vmax=0, cmap="viridis", title="FEE Perfect",
        )
        plt.tight_layout()
        plt.savefig("./FEE.png")
        plt.close()

        fig = plt.subplots(1, 1, figsize=(7, 8))
        bu.plot_healpix(
            data_map=fee_min, vmin=-50, vmax=0, cmap="viridis", title="FEE Minimized",
        )
        plt.tight_layout()
        plt.savefig("./FEE_Min.png")
        plt.close()

        fig = plt.subplots(1, 1, figsize=(7, 8))
        bu.plot_healpix(
            data_map=map_med,
            vmin=-50,
            vmax=0,
            cmap="viridis",
            title=f"{tile} Satellite Map",
        )
        plt.tight_layout()
        plt.savefig(f"./{tile}_Sat.png")
        plt.close()

        fig = plt.subplots(1, 1, figsize=(7, 8))
        data_res = fee - map_med
        fp = fee
        ds = map_med
        er = map_mad
        fp[mask] = np.nan
        ds[mask] = np.nan
        er[mask] = np.nan
        log_chi_sq = np.log(np.nansum(np.square(fp - ds) / er))
        data_res[mask] = np.nan
        bu.plot_healpix(
            data_map=data_res,
            vmin=-5,
            vmax=5,
            cmap="RdYlGn",
            title=f"FEE - {tile} Satellite Map : Log Chisq = {log_chi_sq:.3f}",
        )
        plt.tight_layout()
        plt.savefig(f"./FEE-{tile}_Sat.png")
        plt.close()

        fig = plt.subplots(1, 1, figsize=(7, 8))
        min_res = fee_min - map_med
        fm = fee_min
        ds = map_med
        er = map_mad
        fm[mask] = np.nan
        ds[mask] = np.nan
        er[mask] = np.nan
        log_chi_sq = np.log(np.nansum(np.square(fm - ds) / er))
        min_res[mask] = np.nan
        bu.plot_healpix(
            data_map=min_res,
            vmin=-5,
            vmax=5,
            cmap="RdYlGn",
            title=f"FEE Min - {tile} Satellite Map : Log Chisq = {log_chi_sq:.3f}",
        )
        plt.tight_layout()
        plt.savefig(f"./FEE_Min-{tile}_Sat.png")
        plt.close()
