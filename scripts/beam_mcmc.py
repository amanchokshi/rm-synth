"""Monte Carle to determine MWA beam dipole amplitudes."""

import emcee
import mwa_hyperbeam
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import median_abs_deviation as mad

import beam_utils as bu


def log_likelihood(amps, data):
    """Log likelihood of a beam model give some data.

    Parameter
    ---------
    amps : numpy.array
        16 element dipole amplitude array

    Returns
    -------
    :float
        The log probability that the model fits the data
    """

    # Make a new beam object
    beam = mwa_hyperbeam.FEEBeam()

    # Healpix map with given nside
    nside = 32

    # Zenith angle and Azimuth of healpix pixels
    za, az = bu.healpix_za_az(nside=nside)

    # Satellite beam map frequency
    freq = 138e6

    # Zenith Pointing
    delays = [0] * 16

    # Normalize maps to zenith
    norm_to_zenith = True

    # This is for the trial with 2 parameters
    amps = list(amps) + [1.0] * 14

    # Create model with given amplitudes
    jones = beam.calc_jones_array(az, za, freq, delays, amps, norm_to_zenith)
    unpol_beam = bu.makeUnpolInstrumentalResponse(jones, jones)

    model_XX = np.real(unpol_beam[:, 0])
    #  model_YY = np.real(unpol_beam[:, 3])

    # chisq = np.sum(np.square(data_XX_15 - model_XX) / mad(data_XX_15 - model_XX))
    chisq = np.sum(np.square(data - model_XX))
    log_lik = -0.5 * np.log(chisq)

    return log_lik


def log_prior(amps):
    """Uniform priors in [0, 1]."""

    if np.all(amps <= 1.0) & np.all(amps >= 0.0):
        return 0.0
    else:
        return -np.inf


def log_posterior(amps, data=None):
    """Posterior distribution in log space."""

    # calculate prior
    lp = log_prior(amps)

    # posterior will be -infinity if prior is also -infinity
    if not np.isfinite(lp):
        return -np.inf

    # ln posterior = ln likelihood + ln prior
    return lp + log_likelihood(amps, data)


if __name__ == "__main__":

    # Make a new beam object
    beam = mwa_hyperbeam.FEEBeam()

    # Healpix map with given nside
    nside = 32

    # Zenith angle and Azimuth of healpix pixels
    za, az = bu.healpix_za_az(nside=nside)

    # Satellite beam map frequency
    freq = 138e6

    # Zenith Pointing
    delays = [0] * 16

    # Normalize maps to zenith
    norm_to_zenith = True

    # Create synthetic data at try to recover the input parameters
    amps_15 = [0.7, 0.8] + [1.0] * 14

    jones_15 = beam.calc_jones_array(az, za, freq, delays, amps_15, norm_to_zenith)
    unpol_beam_15 = bu.makeUnpolInstrumentalResponse(jones_15, jones_15)

    data_XX_15 = np.real(unpol_beam_15[:, 0])
    #  data_YY_15 = np.real(unpol_beam_15[:, 3])

    # number of ensemble walkers
    nwalkers = 64

    # Our walkers will be centralised to this location
    #  amps_guess = [0.5] * 16
    amps_guess = [0.5] * 2

    # number of dimensions in sample space
    # ndim = len(amps_guess)
    ndim = 2
    #  ndim = 16

    # Add gaussian perturbation to guess location for each walker
    amps_init = [amps_guess + 1e-1 * np.random.randn(ndim) for i in range(nwalkers)]
    # amps_init = [amps_guess + 1e-1 * np.random.randn(ndim) for i in range(nwalkers)]

    # no. of MCMC iterations - this means there will
    # be n_iterations * nwalkers measurements of the posterior
    n_iterations = 10000

    # Saving MCMC chains
    filename = "beam_mcmc_2amps_v7.h5"
    backend = emcee.backends.HDFBackend(filename)
    backend.reset(nwalkers, ndim)

    # initialise sampler object
    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, log_posterior, kwargs=({"data": data_XX_15}), backend=backend
    )

    # start the chain!
    sampler.run_mcmc(amps_init, n_iterations, progress=True)
