"""Monte Carle to determine MWA beam dipole amplitudes."""

import emcee
import numpy as np

from beam_minimize import beam_mask, likelihood


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

    # number of ensemble walkers
    nwalkers = 64
    ndim = 16

    # Random uniform initial conditions
    amps_init = np.random.rand(nwalkers, ndim)

    # no. of MCMC iterations - this means there will
    # be n_iterations * nwalkers measurements of the posterior
    n_iterations = 200000

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
    )

    # start the chain!
    sampler.run_mcmc(amps_init, n_iterations, progress=True)
