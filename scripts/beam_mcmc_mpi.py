"""Monte Carle to determine MWA beam dipole amplitudes."""

import emcee
import healpy as hp
import mwa_hyperbeam
import numpy as np


def makeUnpolInstrumentalResponse(j1, j2):
    """Convert Jones matracies to unpolarised beam responses.

       Form the visibility matrix in instrumental response from two Jones
       matrices assuming unpolarised sources (hence the brightness matrix is
       the identity matrix)

       Input: j1,j2: Jones matrices of dimension[za][az][2][2]
       Returns: [za][az][[xx,xy],[yx,yy]] where "X" and "Y" are defined by the receptors
       of the Dipole object used in the ApertureArray. Hence to get "XX", you want
       result[za][az][0][0] and for "YY" you want result[za][az][1][1]

       Modified from:
       https://github.com/MWATelescope/mwa_pb/blob/696c2835f44de510da2d5d5bcd3c15223bbe6d7b/mwa_pb/beam_tools.py

       output: [XX, XY, YX, YY]
    """
    result = np.empty_like(j1)

    result[:, 0] = j1[:, 0] * j2[:, 0].conjugate() + j1[:, 1] * j2[:, 1].conjugate()
    result[:, 3] = j1[:, 2] * j2[:, 2].conjugate() + j1[:, 3] * j2[:, 3].conjugate()
    result[:, 1] = j1[:, 0] * j2[:, 2].conjugate() + j1[:, 1] * j2[:, 3].conjugate()
    result[:, 2] = j1[:, 2] * j2[:, 0].conjugate() + j1[:, 3] * j2[:, 1].conjugate()
    return result


def healpix_za_az(nside=32):
    """Zenith angle and Azimuth of healpix map with given nside."""

    npix = hp.nside2npix(nside)

    # healpix indices above horizon
    # convert to zenith angle and azimuth
    above_horizon = range(int(npix / 2))
    za, az = hp.pix2ang(nside, above_horizon)

    return za, az


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
    za, az = healpix_za_az(nside=nside)

    jones_perfect = beam.calc_jones_array(
        az, za, freq, delays, [1.0] * 16, norm_to_zenith
    )
    unpol_perfect = makeUnpolInstrumentalResponse(jones_perfect, jones_perfect)

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
    za, az = healpix_za_az(nside=nside)

    # Indicies of za, az arrays to be masked and not evaluated
    mask_za_az = mask[np.where(mask < az.shape[0])]
    za = np.delete(za, mask_za_az)
    az = np.delete(az, mask_za_az)

    # Create model with given amplitudes
    jones = beam.calc_jones_array(az, za, freq, delays, amps, norm_to_zenith)
    unpol_beam = makeUnpolInstrumentalResponse(jones, jones)

    if pol == "XX":
        model = 10 * np.log10(np.real(unpol_beam[:, 0]))
    else:
        model = 10 * np.log10(np.real(unpol_beam[:, 3]))

    # Remove masked data from sat data maps
    med_map = np.delete(med_map, mask)
    mad_map = np.delete(mad_map, mask)

    chisq = np.sum(np.square(med_map - model) / mad_map)

    return np.log(chisq)


def log_prior(amps):
    """Uniform priors in [0, 1]."""

    if np.all(amps <= 1.0) & np.all(amps >= 0.0):
        return 0.0
    else:
        return -np.inf


def log_posterior(amps, med_map, mad_map, mask, pol):
    """Posterior distribution in log space."""

    # calculate prior
    lp = log_prior(amps)

    # posterior will be -infinity if prior is also -infinity
    if not np.isfinite(lp):
        return -np.inf

    # ln posterior = ln likelihood + ln prior
    return lp + -1 * likelihood(amps, med_map, mad_map, mask, pol)


if __name__ == "__main__":

    import argparse
    from pathlib import Path

    import schwimmbad

    parser = argparse.ArgumentParser(
        description="Determine best fit beam gain parameters - mcmc with mpi.",
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

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--ncores",
        dest="n_cores",
        default=1,
        type=int,
        help="Number of processes (uses multiprocessing).",
    )
    group.add_argument(
        "--mpi", dest="mpi", default=False, action="store_true", help="Run with MPI."
    )

    args = parser.parse_args()

    med_map = np.load(Path(f"{args.map_dir}/rf1_med_maps/{args.map_name}_med.npy"))
    mad_map = np.load(Path(f"{args.map_dir}/rf1_mad_maps/{args.map_name}_mad.npy"))

    if "XX" in args.map_name:
        pol = "XX"
    else:
        pol = "YY"

    mask = beam_mask(med_map, mad_map, pol=pol, db_thresh=-30, zen_mask=20, nside=32)

    pool = schwimmbad.choose_pool(mpi=args.mpi, processes=args.n_cores)

    # number of ensemble walkers
    nwalkers = 64
    ndim = 16

    # Random uniform initial conditions
    amps_init = np.random.rand(nwalkers, ndim)

    # no. of MCMC iterations - this means there will
    # be n_iterations * nwalkers measurements of the posterior
    n_iterations = 20

    # Saving MCMC chains
    out_dir = Path(f"{args.out_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)
    filename = f"{out_dir}/{args.map_name}_{n_iterations}_beam_mcmc.h5"

    backend = emcee.backends.HDFBackend(filename)
    backend.reset(nwalkers, ndim)

    # initialise sampler object
    sampler = emcee.EnsembleSampler(
        nwalkers,
        ndim,
        log_posterior,
        args=(med_map, mad_map, mask, pol),
        backend=backend,
        pool=pool,
    )

    # start the chain!
    sampler.run_mcmc(amps_init, n_iterations, progress=False)
