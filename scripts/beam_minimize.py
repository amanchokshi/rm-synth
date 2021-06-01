"""Minimize beam model to determine MWA beam dipole amplitudes."""

import mwa_hyperbeam
import numpy as np
from scipy.optimize import minimize

import beam_utils as bu


def likelihood(amps, data, mask, pol):
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

    return np.log(chisq)


if __name__ == "__main__":

    import argparse
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

    # Hyperbeam settings
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

    # Our walkers will be centralised to this location
    nwalkers = 512

    # Loop over initial amps and minimize
    min_amps = []

    for i in range(nwalkers):
        print(f"Walker : [{i}/{nwalkers}]")
        result = minimize(
            likelihood,
            np.random.rand(16),
            args=(data_sat, mask_30dB, pol),
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

    out_dir = Path("/astro/mwaeor/achokshi/rm-synth/data/beam_min/")
    out_dir.mkdir(parents=True, exist_ok=True)

    np.save(f"{out_dir}/{map_name}_beam_min_{nwalkers}_walk_mask.npy", min_amps)
