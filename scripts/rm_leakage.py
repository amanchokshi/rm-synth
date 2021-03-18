"""
Plot the 𝜙 = 0 slice of rm cubes as a measure of leakage.
"""

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def read_rm_cube(cube_dir, cube_name):
    """Read rm fits cube created by cuFFS.

    cuFFS creates spectral cubes with the weird
    axis order. The first axis is phi, followed
    by ra, dec. When the fits file is read by
    astropy, the indexing order is reversed

    data == [δ, ⍺, 𝜙]

    Parameters
    ----------
    cube_name : str
        Name of cube fits file
    cube_dir: str
        Path to directory with rm cubes
    """

    # Open cube and grab header and data
    with fits.open(f"{cube_dir}/{cube_name}_p.phi.dirty.fits") as hdus:
        hdu = hdus[0]
        hdr = hdu.header
        data = hdu.data

    # World Coordinate System
    wcs = WCS(hdr)

    # Determine range of phi, ra, dec in data
    # Here, the order of indexing the same as the fits file
    # wcs.all_pix2world(𝜙, ⍺, δ, 0)
    # The last 0 ensures pythonic zero indexing
    phi, _, _ = wcs.all_pix2world(np.arange(data.shape[2]), 0, 0, 0)
    _, ras, _ = wcs.all_pix2world(0, np.arange(data.shape[1]), 0, 0)
    _, _, decs = wcs.all_pix2world(0, 0, np.arange(data.shape[0]), 0)

    phi_0_idx = np.where(phi == 0)[0]

    rm_leak = data[:, :, phi_0_idx]

    plt.style.use("seaborn")
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection=wcs[:, :, int(phi_0_idx)])

    im = ax.imshow(rm_leak, origin="lower", cmap="Spectral_r")

    ax.coords.grid(True, color="white", alpha=0.8, ls="dotted")
    ax.coords[0].set_format_unit(u.deg)
    ax.coords[0].set_auto_axislabel(False)
    ax.set_xlabel("Right Ascension [deg]")

    ax.coords[1].set_format_unit(u.deg)
    ax.coords[1].set_auto_axislabel(False)
    ax.set_ylabel("Declination [deg]")

    #  ax.set_title(f"Quadratic Leakage Surface [{pol}/I]", y=1.2)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="5%", pad=0.1)
    cax.coords[0].set_ticklabel_position("t")
    cax.coords[0].set_axislabel_position("t")

    if "1120300232" in cube_name:
        band = "167-200MHz"
    else:
        band = "200-230MHz"

    if "ana" in cube_name:
        beam = "Analytic"
    else:
        beam = "FEE"

    cax.coords[0].set_axislabel(f"Polarized Flux Density at $\phi$ = 0 : {beam} {band}")
    cax.coords[1].set_ticks(
        alpha=0, color="w", size=0, values=[] * u.dimensionless_unscaled
    )
    plt.colorbar(im, cax=cax, orientation="horizontal")
    plt.show()


if __name__ == "__main__":

    names = ["ana_1120300232", "fee_1120300232", "ana_1120300352", "fee_1120300352"]

    cube_name = names[0]
    cube_dir = "../data/leakage/rm_cubes"

    read_rm_cube(cube_dir, cube_name)
