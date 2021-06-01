"""Minimize beam model to determine MWA beam dipole amplitudes."""

import mwa_hyperbeam
import numpy as np
from scipy.optimize import minimize
#  from matplotlib import pyplot as plt

import beam_utils as bu


def likelihood(amps, data, mask):
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

    # Create model with given amplitudes
    jones = beam.calc_jones_array(az, za, freq, delays, amps, norm_to_zenith)
    unpol_beam = bu.makeUnpolInstrumentalResponse(jones, jones)

    model_XX = 10 * np.log10(np.real(unpol_beam[:, 0]))

    # Remove NaNs from data & model arrays
    model_XX = model_XX[~np.isnan(data)]
    data = data[~np.isnan(data)]

    # Mask nulls
    model_XX = model_XX[mask]
    data = data[mask]

    chisq = np.sum(np.square(data - model_XX))

    return np.log(chisq)


if __name__ == "__main__":

    # Healpix map with given nside
    nside = 32

    # Zenith angle and Azimuth of healpix pixels
    za, az = bu.healpix_za_az(nside=nside)

    # Make a new beam object
    beam = mwa_hyperbeam.FEEBeam()

    # Satellite beam map frequency
    freq = 138e6

    # Zenith Pointing
    delays = [0] * 16

    # Normalize maps to zenith
    norm_to_zenith = True

    data_S06XX = np.load("../data/embers_healpix/S06XX_rf1XX_0.npz")["beam_map"][
        : az.shape[0]
    ]

    jones_perfect = beam.calc_jones_array(
        az, za, freq, delays, [1.0] * 16, norm_to_zenith
    )
    unpol_perfect = bu.makeUnpolInstrumentalResponse(jones_perfect, jones_perfect)
    model_XX_perfect = 10 * np.log10(np.real(unpol_perfect[:, 0]))

    mask_30dB = np.where(model_XX_perfect >= -30)

    # Our walkers will be centralised to this location
    nwalkers = 1024

    # Loop over initial amps and minimize
    min_amps = []

    for i in range(nwalkers):
        print(f"Walker : [{i}/{nwalkers}]")
        result = minimize(
            likelihood,
            np.random.rand(16),
            args=(data_S06XX, mask_30dB),
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

    np.save("S06XX_beam_min_1024_walk_mask.npy", min_amps)
