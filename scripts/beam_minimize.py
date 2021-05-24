"""Minimize beam model to determine MWA beam dipole amplitudes."""

import mwa_hyperbeam
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import minimize
from scipy.stats import median_abs_deviation as mad

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
        The probability that the model fits the data
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

    #  chisq = np.sum(np.square(data - model_XX) / mad(data - model_XX))
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
    #  amps_15 = [0.4] * 1 + [0.8] * 1 + [1.0] * 14
    # amps_15 = [
    #     0.2,
    #     0.3,
    #     0.4,
    #     0.5,
    #     0.6,
    #     0.7,
    #     0.8,
    #     0.9,
    #     0.9,
    #     0.8,
    #     0.7,
    #     0.6,
    #     0.5,
    #     0.4,
    #     0.3,
    #     0.2,
    # ]

    amps_15 = [0.5] + [1.0] * 15

    jones_15 = beam.calc_jones_array(az, za, freq, delays, amps_15, norm_to_zenith)
    unpol_beam_15 = bu.makeUnpolInstrumentalResponse(jones_15, jones_15)

    data_XX_15 = np.real(unpol_beam_15[:, 0])
    # data_YY_15 = np.real(unpol_beam_15[:, 3])

    # lik = []
    # for i in np.linspace(0, 1, 101):
    #     amps = [i] + [1.0] * 15

    #     log_lik = likelihood(amps, data_XX_15)

    #     lik.append(log_lik)

    # plt.style.use("seaborn")
    # plt.plot(np.linspace(0, 1, 101), lik, "-o", color="midnightblue")
    # plt.xlabel("$A_0$ Amplitude")
    # plt.ylabel("Liklihood")
    # plt.tight_layout()
    # plt.show()

    # Our walkers will be centralised to this location
    amps_guess = [0.5] * 16

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
        options={"maxiter": 10000, "disp": True},
    )
    print(result.x)
