"""Monte Carle to determine MWA beam dipole amplitudes."""

import emcee
import mwa_hyperbeam
import numpy as np

import beam_utils as bu


def log_likelihood(amps, data, mask, pol):
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

    # Create model with given amplitudes
    jones = beam.calc_jones_array(az, za, freq, delays, amps, norm_to_zenith)
    unpol_beam = bu.makeUnpolInstrumentalResponse(jones, jones)

    if pol == "XX":
        model = 10 * np.log10(np.real(unpol_beam[:, 0]))
    else:
        model = 10 * np.log10(np.real(unpol_beam[:, 3]))

    # Remove NaNs from data & model arrays
    model = model[~np.isnan(data)]
    data = data[~np.isnan(data)]

    # Mask nulls
    model = model[mask]
    data = data[mask]

    chisq = np.sum(np.square(data - model))

    return -0.5 * np.log(chisq)


def log_prior(amps):
    """Uniform priors in [0, 1]."""

    if np.all(amps <= 1.0) & np.all(amps >= 0.0):
        return 0.0
    else:
        return -np.inf


def log_posterior(amps, data=None, mask=None, pol=None):
    """Posterior distribution in log space."""

    # calculate prior
    lp = log_prior(amps)

    # posterior will be -infinity if prior is also -infinity
    if not np.isfinite(lp):
        return -np.inf

    # ln posterior = ln likelihood + ln prior
    return lp + log_likelihood(amps, data, mask, pol)


if __name__ == "__main__":

    import argparse
    from multiprocessing import Pool
    from pathlib import Path

    parser = argparse.ArgumentParser(
        description="Determine best fit beam gain parameters",
    )

    parser.add_argument(
        "--sat_map",
        metavar="\b",
        type=str,
        required=True,
        help="Path to satellite beam data map",
    )

    args = parser.parse_args()

    sat_map = Path(args.sat_map)

    map_name = sat_map.stem.split("_")[0]

    if "XX" in map_name:
        pol = "XX"
    else:
        pol = "YY"

    # Make a new beam object
    beam = mwa_hyperbeam.FEEBeam()

    # Healpix map with given nside
    nside = 32
    freq = 138e6
    delays = [0] * 16
    norm_to_zenith = True

    # Zenith angle and Azimuth of healpix pixels
    za, az = bu.healpix_za_az(nside=nside)

    # Load satellite beam map
    data_sat = np.load(sat_map)["beam_map"][: az.shape[0]]

    # Create mask based on -30dB threshold of perfect FEE model
    jones_perfect = beam.calc_jones_array(
        az, za, freq, delays, [1.0] * 16, norm_to_zenith
    )
    unpol_perfect = bu.makeUnpolInstrumentalResponse(jones_perfect, jones_perfect)

    if pol == "XX":
        model_perfect = 10 * np.log10(np.real(unpol_perfect[:, 0]))
        mask_30dB = np.where(model_perfect >= -30)
    else:
        model_perfect = 10 * np.log10(np.real(unpol_perfect[:, 3]))
        mask_30dB = np.where(model_perfect >= -30)

    # number of ensemble walkers
    nwalkers = 64
    ndim = 16

    # Add gaussian perturbation to guess location for each walker
    #  amps_init = [amps_guess + 1e-1 * np.random.randn(ndim) for i in range(nwalkers)]
    amps_init = np.random.rand(nwalkers, ndim)

    # no. of MCMC iterations - this means there will
    # be n_iterations * nwalkers measurements of the posterior
    n_iterations = 100000

    # Saving MCMC chains
    filename = f"/astro/mwaeor/achokshi/rm-synth/data/beam_mcmc/beam_mcmc_{map_name}.h5"

    backend = emcee.backends.HDFBackend(filename)
    backend.reset(nwalkers, ndim)

    with Pool() as pool:

        # initialise sampler object
        sampler = emcee.EnsembleSampler(
            nwalkers,
            ndim,
            log_posterior,
            kwargs=({"data": data_sat, "mask": mask_30dB, "pol": pol}),
            backend=backend,
            pool=pool,
        )

        # start the chain!
        sampler.run_mcmc(amps_init, n_iterations, progress=True)
