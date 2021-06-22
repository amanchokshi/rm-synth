"""Evaluate all combinations of gain parameters from beam_minimize.py"""

import concurrent.futures
import itertools
import json
from pathlib import Path

import numpy as np
from scipy import stats
from scipy.signal import find_peaks
from tqdm import tqdm

from beam_minimize import beam_mask, likelihood


def write_json(data, filename=None, out_dir=None):
    """writes data to json file in output dir

    :param data: Data to be written to json file
    :param filename: Json filename :class:`~str`
    :param out_dir: Path to output directory :class:`~str`

    """

    with open(f"{out_dir}/{filename}", "w") as f:
        json.dump(data, f, indent=4)


def peak_amps(beam_min_npy):
    """Find all peaks of diople amplitudes from beam_minimize.py."""

    data = np.load(beam_min_npy)

    amps = []
    for i in range(16):
        kde = stats.gaussian_kde(data[:, i])
        kde_series = kde(np.linspace(0, 1, 2048))
        kde_peak = np.amax(kde_series)
        kde_m = np.append(kde_series, kde_series[-2])
        peaks, _ = find_peaks(kde_m, height=0.3 * kde_peak)
        amps.append(np.linspace(0, 1, 2048)[peaks])

    return amps


def amp_combinations(amps):
    """All possible combinations of peak diople amps."""

    amps_16 = np.array(
        list(
            itertools.product(
                amps[0],
                amps[1],
                amps[2],
                amps[3],
                amps[4],
                amps[5],
                amps[6],
                amps[7],
                amps[8],
                amps[9],
                amps[10],
                amps[11],
                amps[12],
                amps[13],
                amps[14],
                amps[15],
            )
        )
    )

    return amps_16


def amp_comb_chisq(tile):
    """Determine chisq for all combinations of dipole amps."""

    out_dir = Path("../data/beam_min_1024_masked/raw")
    out_dir.mkdir(parents=True, exist_ok=True)

    beam_min = f"../data/beam_min_1024_masked/raw/{tile}_beam_min_1024_walk_mask.npy"

    p_amps = peak_amps(beam_min)
    amps_16 = amp_combinations(p_amps)

    amps_chisq = {}
    for i, a16 in enumerate(tqdm(amps_16)):

        if "XX" in tile:
            pol = "XX"
        else:
            pol = "YY"

        # Load satellite beam map
        map_med = np.load(f"../data/embers_maps/rf1_med_maps/{tile}_med.npy")
        map_mad = np.load(f"../data/embers_maps/rf1_mad_maps/{tile}_mad.npy")

        mask = beam_mask(map_med, map_mad, pol=pol)

        chi = likelihood(a16, map_med, map_mad, mask, pol)

        amps_chisq[i] = [list(a16), chi]

    write_json(amps_chisq, filename=f"{tile}_amp_combinations.json", out_dir=out_dir)
    print(tile)


if __name__ == "__main__":

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

    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(amp_comb_chisq, tiles)

    for result in results:
        print(result)
