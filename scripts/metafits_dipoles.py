"""
Augment metafits files to have non-uniform
dipole gains based on satellite beam maps.
"""

import json

import numpy as np
from astropy.io import fits

np.random.seed(256)


def mod_metafits(meta_dir, obsid, col_name=None, col_data=None, suffix='DipAmps'):
    """Add a DipAmps column to metafits files."""

    # https://docs.astropy.org/en/stable/_modules/astropy/io/fits/column.html
    column = fits.Column(name=col_name, format="16E", array=col_data)
    # column = fits.Column(name=col_name, format="16C", array=col_data)

    with fits.open(f"{meta_dir}/{obsid}.metafits") as hdus:
        table_new = fits.BinTableHDU.from_columns(hdus[1].columns + column)
        hdus[1] = table_new
        hdus.writeto(f"{meta_dir}/{obsid}_{suffix}.metafits")


if __name__ == "__main__":

    # Path to json file it best fit dipole amps from satellite maps
    sat_amps = "../data/beam_min_1024_masked/sat_dipole_amps.json"

    with open(sat_amps, "r") as f:
        data = f.read()
        amps = json.loads(data)

        # Determine XX and YY tile name
        tiles = list(amps.keys())
        XX = [i for i in tiles if "XX" in i]
        YY = [i for i in tiles if "YY" in i]

        # Arrays of XX, YY dipole gains
        amps_xx = np.array([amps[x] for x in XX])
        amps_yy = np.array([amps[y] for y in YY])

    # Random selection of x and y dipole amps for 14 satellite maps
    # Writing out what my past self did -
    # The MCMC found the best parameters for each dipole - 14 tiles, 16 dip, 2 amps
    # Below, for each of the 128 tiles in a full array, for a particular tile,
    # a random value from the 14 possible satellite MCMC values are chosen
    tile_amp_inds = np.random.randint(0, 14, 256).reshape(2, 128)
    mwa_xx_amps = [amps_xx[i] for i in tile_amp_inds[0]]
    mwa_yy_amps = [amps_yy[i] for i in tile_amp_inds[1]]

    # Experimenting with addding 10% phase errors
    # or drawing phase errors from a 0 mean 10% 36deg Normal
    # phase = np.random.normal(
    #     0, scale=2*np.pi/10, size=2*14*16).reshape(2, 14, 16)
    # # phase_xx = phase[0]
    # # phase_yy = phase[1]
    # mwa_xx_amps_phase = [amps_xx[i] *
    #                      np.exp(phase[0][i]*1j) for i in tile_amp_inds[0]]
    # mwa_yy_amps_phase = [amps_yy[i] *
    #                      np.exp(phase[0][i]*1j) for i in tile_amp_inds[1]]

    # The metafits files begin with a row of Y, then X
    mwa_amps = []
    for i in range(128):
        mwa_amps.append(mwa_yy_amps[i])
        mwa_amps.append(mwa_xx_amps[i])

    mwa_amps = np.array(mwa_amps)

    # Do the same for the phase
    # mwa_amps_phase = []
    # for i in range(128):
    #     mwa_amps_phase.append(mwa_yy_amps_phase[i])
    #     mwa_amps_phase.append(mwa_xx_amps_phase[i])
    #
    # mwa_amps_phase = np.array(mwa_amps_phase)

    # Modify and write new metafits file

    # obsids = [1120082744, 1120082864]
    obsids = [1088285600]
    # meta_dir = "../data/metafits_dipoles"
    meta_dir = "."

    # print(f"Saving Modified Metafits to : {meta_dir}")
    # for obs in obsids:
    #     mod_metafits(meta_dir, obs, col_name="DipAmps", col_data=mwa_amps)
    #     # mod_metafits(meta_dir, obs, col_name="DipAmps", col_data=mwa_amps_phase)
    #     print(f"Augmented {obs}.metafits with dipole amps")

    # In this section we only change the dipole amps of one tile
    # Tile036
    with fits.open(f'{obsids[0]}.metafits') as hdul:
        hdr = hdul[1].header
        data = hdul[1].data

        # Find indicies of dipoles in Tile036
        T36_idxs = np.where(data['TileName'] == 'Tile036')

        # Create a dipole array with only Tile036 deformed
        T36 = np.zeros(mwa_amps.shape)
        T36[T36_idxs] += mwa_amps[T36_idxs]
        T36 = np.where(T36 == 0, 1, T36)

    print(f"Saving Modified Metafits to : {meta_dir}")
    for obs in obsids:
        mod_metafits(meta_dir, obs, col_name="DipAmps", col_data=T36, suffix='Tile036')
        print(f"Augmented {obs}.metafits with dipole amps")
