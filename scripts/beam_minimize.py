"""Minimize beam model to determine MWA beam dipole amplitudes."""

import healpy as hp
import mwa_hyperbeam
import numpy as np
from scipy.optimize import minimize

import beam_utils as bu


def beam_mask(med_map, mad_map, pol=None, db_thresh=-30, zen_mask=20, nside=32):
    """Create the ultimate beam mask.

    The satellite beam map needs to be maskes in various ways before
    it can be used for minimization.

     - Mask the nulls below `db_thresh` from zenith power
     - Mask the central `zen_mask` degrees
     - Mask NaNs in the Median satellite map
     - Mask zeros in MAD map - indicative of single satellite pass in the pixel

    Parameter
    ---------
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
    beam = mwa_hyperbeam.FEEBeam()

    # Hyperbeam settings
    freq = 138e6
    delays = [0] * 16
    norm_to_zenith = True

    # Zenith angle and Azimuth of healpix pixels
    za, az = bu.healpix_za_az(nside=nside)

    jones_perfect = beam.calc_jones_array(
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


def likelihood(amps, med_map, mad_map, mask, pol):
    """Likelihood of a beam model give some data.

    Parameter
    ---------
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
    beam = mwa_hyperbeam.FEEBeam()

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
    jones = beam.calc_jones_array(az, za, freq, delays, amps, norm_to_zenith)
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


if __name__ == "__main__":

    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(
        description="Determine best fit beam gain parameters",
    )

    parser.add_argument(
        "--map_dir",
        metavar="\b",
        type=str,
        required=True,
        help="Path to satellite map directory",
    )

    parser.add_argument(
        "--map_name",
        metavar="\b",
        type=str,
        required=True,
        help="Name of satellite map - Ex: S06XX_rf1XX",
    )

    parser.add_argument(
        "--out_dir",
        metavar="\b",
        type=str,
        required=True,
        help="Output directory to save minimization results",
    )

    args = parser.parse_args()

    med_map = np.load(Path(f"{args.map_dir}/rf1_med_maps/{args.map_name}_med.npy"))
    mad_map = np.load(Path(f"{args.map_dir}/rf1_mad_maps/{args.map_name}_mad.npy"))

    if "XX" in args.map_name:
        pol = "XX"
    else:
        pol = "YY"

    mask = beam_mask(med_map, mad_map, pol=pol, db_thresh=-30, zen_mask=20, nside=32)

    # Our walkers will be centralised to this location
    nwalkers = 1024

    # Loop over initial amps and minimize
    min_amps = []

    for i in range(nwalkers):
        print(f"Walker : [{i}/{nwalkers}]")
        result = minimize(
            likelihood,
            np.random.rand(16),
            args=(med_map, mad_map, mask, pol),
            bounds=(
                (0, 1),
                (0, 1),
                (0, 1),
                (0, 1),
                (0, 1),
                (0, 1),
                (0, 1),
                (0, 1),
                (0, 1),
                (0, 1),
                (0, 1),
                (0, 1),
                (0, 1),
                (0, 1),
                (0, 1),
                (0, 1),
            ),
            options={"maxiter": 10000, "disp": False},
        )
        min_amps.append(result.x)
        #  print(result.x)

    min_amps = np.array(min_amps)

    out_dir = Path(f"{args.out_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)

    np.save(f"{out_dir}/{args.map_name}_beam_min_{nwalkers}_walk_mask.npy", min_amps)
