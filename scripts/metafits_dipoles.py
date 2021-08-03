"""
Augment metafits files to have non-uniform
dipole gains based on satellite beam maps.
"""

import json

import numpy as np
from astropy.io import fits

np.random.seed(256)


def mod_metafits(meta_dir, obsid, col_name=None, col_data=None):
    """Add a DipAmps column to metafits files."""

    # https://docs.astropy.org/en/stable/_modules/astropy/io/fits/column.html
    column = fits.Column(name=col_name, format="16E", array=col_data)

    with fits.open(f"{meta_dir}/{obsid}.metafits") as hdus:
        table_new = fits.BinTableHDU.from_columns(hdus[1].columns + column)
        hdus[1] = table_new
        hdus.writeto(f"{meta_dir}/{obsid}_DipAmps.metafits")


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
    tile_amp_inds = np.random.randint(0, 14, 256).reshape(2, 128)
    mwa_xx_amps = [amps_xx[i] for i in tile_amp_inds[0]]
    mwa_yy_amps = [amps_yy[i] for i in tile_amp_inds[1]]

    # The metafits files begin with a row of Y, then X
    mwa_amps = []
    for i in range(128):
        mwa_amps.append(mwa_yy_amps[i])
        mwa_amps.append(mwa_xx_amps[i])

    mwa_amps = np.array(mwa_amps)

    # Modify and write new metafits file

    obsids = [1120082744, 1120082864]
    meta_dir = "../data/metafits_dipoles"

    print(f"Saving Modified Metafits to : {meta_dir}")
    for obs in obsids:
        mod_metafits(meta_dir, obs, col_name="DipAmps", col_data=mwa_amps)
        print(f"Augmented {obs}.metafits with dipole amps")
