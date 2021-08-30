import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import pyplot as plt


def read_rm_cube(rm_coords, rm_fits_cube):
    """Read rm fits cube created by cuFFS.

    cuFFS creates spectral cubes with the weird
    axis order. The first axis is phi, followed
    by ra, dec. When the fits file is read by
    astropy, the indexing order is reversed

    data == [δ, ⍺, 𝜙]

    Parameters
    ----------
    rm_coords : skycoord object
        SkyCoord object of RM source
    rm_fits_cube: str
        Path to rm fits cube
    """

    # Open cube and grab header and data
    with fits.open(rm_fits_cube) as hdus:
        hdu = hdus[0]
        hdr = hdu.header
        data_p = hdu.data

    # World Coordinate System
    wcs = WCS(hdr)

    # Determine range of phi, ra, dec in data
    # Here, the order of indexing the same as the fits file
    # wcs.all_pix2world(𝜙, ⍺, δ, 0)
    # The last 0 ensures pythonic zero indexing
    phi, _, _ = wcs.all_pix2world(np.arange(data_p.shape[2]), 0, 0, 0)
    _, ras, _ = wcs.all_pix2world(0, np.arange(data_p.shape[1]), 0, 0)
    _, _, decs = wcs.all_pix2world(0, 0, np.arange(data_p.shape[0]), 0)

    # Determine pixel coordinates of pogs object
    _, ra_x, dec_y = wcs.all_world2pix(0, rm_coords.ra.deg, rm_coords.dec.deg, 0)

    # Convert to closes integer pixel index
    ra_x = np.round(ra_x).astype(int)
    dec_y = np.round(dec_y).astype(int)

    # Extract spectral slice of cube at location of pogs source
    spec = data_p[np.round(dec_y).astype(int), np.round(ra_x).astype(int), :]

    # Determine RM peak and index in cube
    spec_peak = np.amax(spec)
    phi_z = np.where(spec == spec_peak)[0][0]

    return ra_x, dec_y, phi_z, phi, data_p, wcs


if __name__ == "__main__":

    # Epic spectral colours
    colours = ["#DA3752",
               "#FCAD61",
               "#66C1A4",
               "#3287BC",
               "#5E4FA1"]

    rm_coords = SkyCoord(
        ra=10.44137913863797 * u.degree, dec=-26.78792973842179 * u.degree, frame="icrs"
    )
    rm_cube = "../data/1120082744/1120082744_dip_0/imgs/cubes/1120082744_dip_0_1120082744_p.phi.dirty.fits"

    ra_x, dec_y, phi_z, phi, data_p, wcs = read_rm_cube(rm_coords, rm_cube)

    rm_spec = data_p[dec_y, ra_x, :]
    arr = [
        data_p[dec_y, ra_x, :],
        data_p[dec_y + 1, ra_x, :],
        data_p[dec_y + 1, ra_x + 1, :],
        data_p[dec_y + 1, ra_x - 1, :],
        data_p[dec_y, ra_x + 1, :],
        data_p[dec_y, ra_x - 1, :],
        data_p[dec_y - 1, ra_x, :],
        data_p[dec_y - 1, ra_x + 1, :],
        data_p[dec_y - 1, ra_x - 1, :],
    ]

    rm_med = np.median(arr, axis=0)
    rm_max = np.amax(arr, axis=0)

    plt.rcParams.update(
        {
            #  "font.size": 15,
            #  "text.usetex": True,
            "font.family": "serif",
            #  "font.serif": "Times New Roman",
            # Use 10pt font in plots, to match 10pt font in document
            "axes.labelsize": 8,
            "axes.titlesize": 9,
            "font.size": 8,
            # Make the legend/label fonts a little smaller
            "legend.fontsize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
        }
    )
    plt.style.use("seaborn")
    plt.plot(phi, rm_max, color=colours[0], linewidth=2, alpha=0.9)
    plt.xlim([-200, 200])
    plt.xlabel("Faraday Depth [rad/m$^2$]")
    plt.ylabel("Polarized Flux Density [Jy/PSF/RMSF]")
    plt.tight_layout()
    plt.savefig("../data/1120082744/1120082744_dip_0/imgs/cubes/RM_1120082744_dip_0.pdf")
    #  plt.show()
