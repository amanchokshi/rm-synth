"""Evaluate all combinations of gain parameters from beam_minimize.py"""

import concurrent.futures
import itertools
import json
from pathlib import Path

import healpy as hp
import mwa_hyperbeam
import numpy as np
from scipy import stats
from scipy.signal import find_peaks
from tqdm import tqdm

import beam_utils as bu


def write_json(data, filename=None, out_dir=None):
    """writes data to json file in output dir

    :param data: Data to be written to json file
    :param filename: Json filename :class:`~str`
    :param out_dir: Path to output directory :class:`~str`

    """

    with open(f"{out_dir}/{filename}", "w") as f:
        json.dump(data, f, indent=4)


def beam_mask(
    hyperbeam, med_map, mad_map, pol=None, db_thresh=-30, zen_mask=20, nside=32
):
    """Create the ultimate beam mask.

    The satellite beam map needs to be maskes in various ways before
    it can be used for minimization.

     - Mask the nulls below `db_thresh` from zenith power
     - Mask the central `zen_mask` degrees
     - Mask NaNs in the Median satellite map
     - Mask zeros in MAD map - indicative of single satellite pass in the pixel

    Parameter
    ---------
    hyperbeam : :class:`mwa_hyperbeam.FEEBeam()`
        Instance of FEEBeam class
    med_map : numpy.array
        Healpix satellite median map
    mad_map : numpy.array
        Healpix satellite MAD map - noise
    pol : string
        Polarization of map - XX / YY
    db_thresh : int / float
        Null mask power threshold - default: -30dB
    zen_mask : int / float
        Central radii in degrees to mask - default: 20
    nside : int
        Healpix nside - default: 32

    Returns
    -------
    :numpy.array
        Numpy array of indicies of nside healpix map to be masked

    """

    # Make a new beam object
    # beam = mwa_hyperbeam.FEEBeam()

    # Hyperbeam settings
    freq = 138e6
    delays = [0] * 16
    norm_to_zenith = True

    # Zenith angle and Azimuth of healpix pixels
    za, az = bu.healpix_za_az(nside=nside)

    jones_perfect = hyperbeam.calc_jones_array(
        az, za, freq, delays, [1.0] * 16, norm_to_zenith
    )
    unpol_perfect = bu.makeUnpolInstrumentalResponse(jones_perfect, jones_perfect)

    if pol == "XX":
        fee_hyperbeam = 10 * np.log10(np.real(unpol_perfect[:, 0]))
    else:
        fee_hyperbeam = 10 * np.log10(np.real(unpol_perfect[:, 3]))

    # Anything below 30dB from peak of FEE zenith norm beam
    mask_dB = np.where(fee_hyperbeam < db_thresh)[0]

    # The central pixels upto `zen_mask` degrees
    zenith_mask = np.arange(hp.ang2pix(nside, np.deg2rad(zen_mask), 0))

    # All nans in sat median map
    nan_mask = np.where(np.isnan(med_map) == True)[0]

    # All 0 values in MAD array - indicative of only single satellite pass
    mad_mask = np.where(mad_map == 0.0)[0]

    mask = np.unique(np.concatenate((zenith_mask, mask_dB, nan_mask, mad_mask)))

    return mask


def likelihood(hyperbeam, amps, med_map, mad_map, mask, pol):
    """Likelihood of a beam model give some data.

    Parameter
    ---------
    hyperbeam : :class:`mwa_hyperbeam.FEEBeam()`
        Instance of FEEBeam class
    amps : numpy.array
        16 element dipole amplitude array
    data : numpy.array
        Satellite beam model healpix array
    mask : numpy.array
        Indicies of data array where beam power >= -30dB
    pol : string
        Either XX, YY, indicating the polarization of current map

    Returns
    -------
    :float
        The log probability that the model fits the data
    """

    # Make a new beam object
    # beam = mwa_hyperbeam.FEEBeam()

    # Hyperbeam settings
    nside = 32
    freq = 138e6
    delays = [0] * 16
    norm_to_zenith = True

    # Zenith angle and Azimuth of healpix pixels
    za, az = bu.healpix_za_az(nside=nside)

    # Indicies of za, az arrays to be masked and not evaluated
    mask_za_az = mask[np.where(mask < az.shape[0])]
    za = np.delete(za, mask_za_az)
    az = np.delete(az, mask_za_az)

    # Create model with given amplitudes
    jones = hyperbeam.calc_jones_array(az, za, freq, delays, amps, norm_to_zenith)
    unpol_beam = bu.makeUnpolInstrumentalResponse(jones, jones)

    if pol == "XX":
        model = 10 * np.log10(np.real(unpol_beam[:, 0]))
    else:
        model = 10 * np.log10(np.real(unpol_beam[:, 3]))

    # Remove masked data from sat data maps
    med_map = np.delete(med_map, mask)
    mad_map = np.delete(mad_map, mask)

    chisq = np.sum(np.square(med_map - model) / mad_map)

    return np.log(chisq)


def peak_amps(beam_min_npy):
    """Find all peaks of diople amplitudes from beam_minimize.py."""

    data = np.load(beam_min_npy)

    amps = []
    for i in range(16):

        if np.unique(data[:, i]).shape[0] == 1:
            amps.append(np.unique(data[:, i]))
        else:
            kde = stats.gaussian_kde(data[:, i], bw_method=0.1)
            kde_series = kde(np.linspace(0, 1, 2048))
            kde_peak = np.amax(kde_series)
            kde_m = np.append(kde_series, kde_series[-2])
            peaks, _ = find_peaks(kde_m, height=0.5 * kde_peak)

            # Sort peaks and pick top two
            peak_height_sort = np.array(sorted(zip(kde_m[peaks], peaks), reverse=True))

            if peak_height_sort.shape[0] > 1:
                amps.append(
                    np.linspace(0, 1, 2048)[peak_height_sort[:2, 1].astype(int)]
                )
            else:
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

    print(tile)
    out_dir = Path("../data/beam_min_1024_masked/raw")
    out_dir.mkdir(parents=True, exist_ok=True)

    beam_min = f"../data/beam_min_1024_masked/raw/{tile}_beam_min_1024_walk_mask.npy"

    # Make a new beam object
    hyperbeam = mwa_hyperbeam.FEEBeam()

    p_amps = peak_amps(beam_min)
    print(p_amps)
    amps_16 = amp_combinations(p_amps)

    amps_chisq = {}
    for i, a16 in enumerate(tqdm(amps_16)):
        # for i, a16 in enumerate(amps_16):

        if "XX" in tile:
            pol = "XX"
        else:
            pol = "YY"

        # Load satellite beam map
        map_med = np.load(f"../data/embers_maps/rf1_med_maps/{tile}_med.npy")
        map_mad = np.load(f"../data/embers_maps/rf1_mad_maps/{tile}_mad.npy")

        mask = beam_mask(hyperbeam, map_med, map_mad, pol=pol)

        chi = likelihood(hyperbeam, a16, map_med, map_mad, mask, pol)

        amps_chisq[i] = [list(a16), chi]

    del hyperbeam

    write_json(amps_chisq, filename=f"{tile}_amp_combinations.json", out_dir=out_dir)


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Determine best fit beam gain parameters",
    )

    parser.add_argument(
        "--tile",
        metavar="\b",
        type=str,
        required=True,
        help="Tile name. Ex: S06XX_rf1XX",
    )

    args = parser.parse_args()

    # tiles = [
    #     "S06XX_rf1XX",
    #     "S06YY_rf1YY",
    #     "S07XX_rf1XX",
    #     "S07YY_rf1YY",
    #     "S08XX_rf1XX",
    #     "S08YY_rf1YY",
    #     "S09XX_rf1XX",
    #     "S09YY_rf1YY",
    #     "S10XX_rf1XX",
    #     "S10YY_rf1YY",
    #     "S12XX_rf1XX",
    #     "S12YY_rf1YY",
    #     "S29XX_rf1XX",
    #     "S29YY_rf1YY",
    #     "S30XX_rf1XX",
    #     "S30YY_rf1YY",
    #     "S31XX_rf1XX",
    #     "S31YY_rf1YY",
    #     "S32XX_rf1XX",
    #     "S32YY_rf1YY",
    #     "S33XX_rf1XX",
    #     "S33YY_rf1YY",
    #     "S34XX_rf1XX",
    #     "S34YY_rf1YY",
    #     "S35XX_rf1XX",
    #     "S35YY_rf1YY",
    #     "S36XX_rf1XX",
    #     "S36YY_rf1YY",
    # ]

    # with concurrent.futures.ProcessPoolExecutor(max_workers=7) as executor:
    # with concurrent.futures.ProcessPoolExecutor() as executor:
    #     results = executor.map(amp_comb_chisq, tiles)

    try:
        amp_comb_chisq(args.tile)
    except Exception as e:
        print(e)
