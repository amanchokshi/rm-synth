"""Minimize beam model to determine MWA beam dipole amplitudes."""

import mwa_hyperbeam
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import minimize

import beam_utils as bu


def likelihood(amps, data):
    """Likelihood of a beam model give some data.

    Parameter
    ---------
    amps : numpy.array
        16 element dipole amplitude array

    Returns
    -------
    :float
        The inverse log probability that the model fits the data
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
    #  amps = list(amps) + [1.0] * 14

    # Create model with given amplitudes
    jones = beam.calc_jones_array(az, za, freq, delays, amps, norm_to_zenith)
    unpol_beam = bu.makeUnpolInstrumentalResponse(jones, jones)

    model_XX = np.real(unpol_beam[:, 0])
    #  model_YY = np.real(unpol_beam[:, 3])

    chisq = np.sum(np.square(data - model_XX))

    return np.log(chisq)


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
    amps_15 = [
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.8,
        0.9,
        0.9,
        0.8,
        0.7,
        0.6,
        0.5,
        0.4,
        0.3,
        0.2,
    ]

    jones_15 = beam.calc_jones_array(az, za, freq, delays, amps_15, norm_to_zenith)
    unpol_beam_15 = bu.makeUnpolInstrumentalResponse(jones_15, jones_15)

    data_XX_15 = np.real(unpol_beam_15[:, 0])
    # data_YY_15 = np.real(unpol_beam_15[:, 3])

    # Our walkers will be centralised to this location
    nwalkers = 64
    amps_guess = [0.5] * 16
    amps_init = [
        amps_guess + 1e-1 * np.random.randn(len(amps_guess)) for i in range(nwalkers)
    ]

    # Loop over initial amps and minimize

    min_amps = []

    for i in range(nwalkers):
        result = minimize(
            likelihood,
            amps_guess,
            args=(data_XX_15),
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

    min_amps = np.array(min_amps)

    np.save("beam_min.npy", min_amps)
