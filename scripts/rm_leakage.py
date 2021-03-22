"""
Plot the 𝜙 = 0 slice of rm cubes as a measure of leakage.
"""

from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

try:
    from skyfield.api import wgs84

    mwa_loc = wgs84.latlon(-26.703319, 116.670815, 337.83)

except Exception as e:
    print(e)
    from skyfield.api import Topos

    mwa_loc = Topos(latitude=-26.703319, longitude=116.670815, elevation_m=337.83)


from instru_leakage import get_beam_weights, read_metafits


def read_rm_cube(rm_cube, noise_cube, metafits, obsid, beam, band, fname, outdir):
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

    # Read noise image
    with fits.open(noise_cube) as hdus:
        noise = hdus[0].data

    # Open cube and grab header and data
    with fits.open(rm_cube) as hdus:
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

    # Pixel indices along Ra: x, Dec: y
    x = np.arange(hdr["NAXIS2"])
    y = np.arange(hdr["NAXIS3"])

    # Pixel meshgrid
    X, Y = np.meshgrid(x, y)
    _, ras, decs, = wcs.wcs_pix2world(phi_0_idx, X, Y, 0)

    # Rescale ras from -180, 180 to 0, 360
    #  ras = np.mod(ras, 360)

    # Flatten 2D ras, decs arrays
    ra_f = ras.flatten()
    dec_f = decs.flatten()

    lst, obsid, freqcent, ra_point, dec_point, delays = read_metafits(metafits)
    beam_weights = get_beam_weights(
        ras=ra_f, decs=dec_f, LST=lst, mwa_lat=mwa_lat, freqcent=freqcent, delays=delays
    )
    beam_weights = beam_weights.reshape(ras.shape)

    # Plotting stuff
    plt.style.use("seaborn")
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection=wcs[:, :, int(phi_0_idx)])

    im = ax.imshow(
        rm_leak[:, :, 0] - noise, origin="lower", cmap="Spectral_r", vmin=0.0, vmax=0.12
    )

    levels = [0.001, 0.01, 0.1, 0.3, 0.6, 0.9]
    CS = ax.contour(
        beam_weights.real,
        levels,
        colors="#222222",
        linewidths=2.1,
        linestyles="dotted",
    )
    ax.clabel(CS, inline=1, fontsize=9)

    ax.coords.grid(True, color="white", alpha=1, linewidth=2.1, ls="dotted")
    ax.coords[0].set_format_unit(u.deg)
    ax.coords[0].set_auto_axislabel(False)
    ax.set_xlabel("Right Ascension [deg]")

    ax.coords[1].set_format_unit(u.deg)
    ax.coords[1].set_auto_axislabel(False)
    ax.set_ylabel("Declination [deg]")

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="5%", pad=0.1)
    cax.coords[0].set_ticklabel_position("t")
    cax.coords[0].set_axislabel_position("t")

    cax.coords[0].set_axislabel(
        f"Polarized Flux at $\phi$ = 0: {obsid}: {beam.upper()}: {band}: [{ra_point:.1f} {dec_point:.1f}]"
    )
    cax.coords[1].set_ticks(
        alpha=0, color="w", size=0, values=[] * u.dimensionless_unscaled
    )
    plt.colorbar(im, cax=cax, orientation="horizontal")
    plt.savefig(f"{outdir}/{fname}", bbox_inches="tight", dpi=300)


if __name__ == "__main__":

    # MWA coordinates
    mwa_lat = mwa_loc.latitude.degrees
    mwa_lon = mwa_loc.longitude.degrees
    mwa_el = mwa_loc.elevation.m

    low_band = ["1120300232", "1120082744"]
    high_band = ["1120300352", "1120082864"]

    obsids = low_band + high_band

    for o in obsids:

        for b in ["ana", "fee"]:

            rm_cube = f"../data/{o}/{b}_wide/imgs/cubes/{b}_wide_{o}_p.phi.dirty.fits"
            noise_cube = f"../data/{o}/{b}_wide/imgs/cubes/{b}_wide_{o}_cube_noise.fits"
            metafits = Path(f"../data/{o}/{b}_wide/{o}_metafits_ppds.fits")

            print(f" ** INFO: Crunching data - {rm_cube}")

            if o in low_band:
                band = "167-200 MHz"
                fname = f"{o}_{b.upper()}_167-200MHz_rm_phi_0.png"
            else:
                band = "200-230 MHz"
                fname = f"{o}_{b.upper()}_200-230MHz_rm_phi_0.png"

            read_rm_cube(rm_cube, noise_cube, metafits, o, b, band, fname, "./")
