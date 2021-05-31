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
    #  norm_to_zenith = False

    # Create model with given amplitudes
    jones = beam.calc_jones_array(az, za, freq, delays, amps, norm_to_zenith)
    unpol_beam = bu.makeUnpolInstrumentalResponse(jones, jones)

    model_XX = 10 * np.log10(np.real(unpol_beam[:, 0]))
    #  model_YY = 10 * np.log10(np.real(unpol_beam[:, 3]))

    # Remove NaNs from data & model arrays
    model_XX = model_XX[~np.isnan(data)]
    #  model_YY = model_YY[~np.isnan(data)]
    data = data[~np.isnan(data)]

    chisq = np.sum(np.square(data - model_XX))
    #  chisq = np.sum(np.square(data - model_YY))

    return np.log(chisq)


if __name__ == "__main__":

    # Healpix map with given nside
    nside = 32

    # Zenith angle and Azimuth of healpix pixels
    za, az = bu.healpix_za_az(nside=nside)

    ###################################################################

    # Make a new beam object
    beam = mwa_hyperbeam.FEEBeam()

    # Satellite beam map frequency
    freq = 138e6

    # Zenith Pointing
    delays = [0] * 16

    # Normalize maps to zenith
    norm_to_zenith = True

    # Create synthetic data at try to recover the input parameters
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
    # amps_15 = np.linspace(0.85, 1.0, 16)
    #  amps_15 = [1.0, 0.95] + [1.0] * 13 + [0.97]

    #  jones_15 = beam.calc_jones_array(az, za, freq, delays, amps_15, norm_to_zenith)
    #  unpol_beam_15 = bu.makeUnpolInstrumentalResponse(jones_15, jones_15)

    #  data_XX_15 = 10 * np.log10(np.real(unpol_beam_15[:, 0]))
    # data_YY_15 = np.real(unpol_beam_15[:, 3])

    data_S06XX = np.load("../data/embers_healpix/S06XX_rf1XX_0.npz")["beam_map"][
        : az.shape[0]
    ]
    data_S06YY = np.load("../data/embers_healpix/S06YY_rf1YY_0.npz")["beam_map"][
        : az.shape[0]
    ]

    #  chi = []
    #  for i in np.linspace(0.0, 1.0, 101):
    #  amps = [i] * 16
    #  prob = likelihood(amps, data_S06XX)
    #  chi.append(prob)

    #  plt.style.use("seaborn")
    #  plt.plot(np.linspace(0.0, 1.0, 101), chi, "-o", color="midnightblue")
    #  plt.xlabel("$A_0$")
    #  plt.ylabel("Prob")
    #  plt.tight_layout()
    #  plt.show()

    ###################################################################

    # Our walkers will be centralised to this location
    nwalkers = 1024
    # amps_guess = [0.5] * 16
    # amps_init = [
    #     amps_guess + 1e-1 * np.random.randn(len(amps_guess)) for i in range(nwalkers)
    # ]

    # Loop over initial amps and minimize
    min_amps = []

    for i in range(nwalkers):
        print(f"Walker : [{i}/{nwalkers}]")
        result = minimize(
            likelihood,
            np.random.rand(16),
            args=(data_S06XX),
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
            #  method="Powell",
        )
        min_amps.append(result.x)
        #  print(result.x)

    min_amps = np.array(min_amps)

    np.save("S06XX_beam_min_1024_walk.npy", min_amps)
