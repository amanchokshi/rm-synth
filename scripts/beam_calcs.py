import mwa_hyperbeam
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm


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

    # Create a grid across the sky at which to evaluate the beam
    nside = 1001

    x_range = np.linspace(-np.pi / 2, np.pi / 2, nside)
    y_range = np.linspace(-np.pi / 2, np.pi / 2, nside)

    x_mesh, y_mesh = np.meshgrid(x_range, y_range)

    az_grid = np.pi - np.arctan2(x_mesh, y_mesh)
    za_grid = np.sqrt(x_mesh ** 2 + y_mesh ** 2)

    below_horizon = np.where(np.sqrt(x_mesh ** 2 + y_mesh ** 2) > np.pi / 2)

    az_grid[below_horizon] = np.nan
    za_grid[below_horizon] = np.nan

    az_f = az_grid.flatten()
    za_f = za_grid.flatten()

    freq = 167000000
    delays = [0] * 16
    amps = [1.0] * 16
    norm_to_zenith = True

    # beam.calc_jones is also available, but that makes a single Jones matrix at a
    # time, so one would need to iterate over az and za. calc_jones_array is done in
    # parallel with Rust (so it's fast).
    jones = beam.calc_jones_array(az_f, za_f, freq, delays, amps, norm_to_zenith)

    unpol_beam = makeUnpolInstrumentalResponse(jones, jones)

    XX = unpol_beam[:, 0]
    YY = unpol_beam[:, 3]

    beam_weights = (XX + YY) / 2.0

    beam_weights = beam_weights.reshape(az_grid.shape)

    #  plt.style.use("seaborn")
    #  plt.imshow(beam_weights.real, norm=LogNorm(vmin=0.0001, vmax=1))
    #  plt.tight_layout()
    #  plt.show()

    levels = [0.0001, 0.001, 0.01, 0.1, 0.3, 0.9]

    fig, ax = plt.subplots()
    im = ax.imshow(beam_weights.real, norm=LogNorm(vmin=0.0001, vmax=1))
    CS = ax.contour(beam_weights.real, levels)
    ax.clabel(CS, inline=1, fontsize=6)
    plt.tight_layout()
    plt.show()
