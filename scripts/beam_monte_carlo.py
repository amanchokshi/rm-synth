"""Monte Carle to determine MWA beam dipole amplitudes."""

import healpy as hp
import mwa_hyperbeam
import numpy as np
from matplotlib import pyplot as plt
from scipy import random
from tqdm import tqdm

import healpix as hpx


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


if __name__ == "__main__":

    # We can make a new beam object with a path to the HDF5 file specified by MWA_BEAM_FILE.
    beam = mwa_hyperbeam.FEEBeam()

    # Healpix map with given nside
    nside = 32
    za, az = hpx.healpix_za_az(nside=nside)

    freq = 138e6
    delays = [0] * 16
    norm_to_zenith = True

    # Setup a beam with the first dipole dead
    # We'll try and determine which dipole is dead using monte carlo magic
    amps_15 = [0.0] * 1 + [1.0] * 15
    jones_15 = beam.calc_jones_array(az, za, freq, delays, amps_15, norm_to_zenith)
    unpol_beam_15 = makeUnpolInstrumentalResponse(jones_15, jones_15)

    XX_15 = np.real(unpol_beam_15[:, 0])
    YY_15 = np.real(unpol_beam_15[:, 3])

    # Monte Carlo Magic Below

    # define a set of N random dipole amps between 0, 1
    N = 1000
    amps_rand = random.uniform(low=0.0, high=1.0, size=(16, N))

    for i in tqdm(range(N)):
        amps = amps_rand[:, i]

        jones = beam.calc_jones_array(az, za, freq, delays, amps, norm_to_zenith)
        unpol_beam = makeUnpolInstrumentalResponse(jones, jones)

        XX = np.real(unpol_beam[:, 0])
        YY = np.real(unpol_beam[:, 3])

    #  npix = hp.nside2npix(nside=32)
    #  beam_response_XX = np.zeros(npix)
    #  beam_response_YY = np.zeros(npix)
    #  above_horizon = range(int(npix / 2))
    #  beam_response_XX[above_horizon] = XX
    #  beam_response_YY[above_horizon] = YY

    #  plt.style.use("seaborn")
    #  fig = plt.figure(figsize=(6, 6))
    #  fig.suptitle("MWA HYPERBEAM FEE MAP XX", fontsize=16, y=0.92)
    #  hpx.plot_healpix(
    #  data_map=10 * np.log10(beam_response_XX),
    #  sub=(1, 1, 1),
    #  cmap="viridis",
    #  vmin=-50,
    #  vmax=0,
    #  )
    #  plt.show()
    #  #  plt.savefig("XX.png", bbox_inches="tight")
