"""
Plot the RM Specta and RMSF of a POGS source.
"""

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from matplotlib import pyplot as plt


def pogs_obj_loc(pogs_obj, pogs_path):
    """Read the POGs tables and create skycoord object for a source

    Parameters
    ----------
    pogs_obj : str
        Name of object in POGS catalogue
    pogs_path : str
        Path to POGS catalogue fits file

    Returns
    -------
    SkyCoord object
        Astropy.coordinates.SkyCoord object
    """

    exgal = Table.read(pogs_path)
    df_eg = exgal.to_pandas()
    df_eg["catalog_id"] = df_eg["catalog_id"].str.decode("utf-8")
    df_eg.set_index("catalog_id", inplace=True)

    pogs_source = df_eg.loc[f"{pogs_obj}"]
    pogs_pos = SkyCoord(pogs_source.ra, pogs_source.dec, unit="deg")

    return pogs_pos


def read_rm_cube(pogs_pos, rm_prefix, cube_dir):
    """Read rm fits cube created by cuFFS.

    cuFFS creates spectral cubes with the weird
    axis order. The first axis is phi, followed
    by ra, dec. When the fits file is read by
    astropy, the indexing order is reversed

    data == [δ, ⍺, 𝜙]

    Parameters
    ----------
    pogs_pos : skycoord object
        SkyCoord object from pogs_obj_loc
    rm_prefix : str
        Prefix to cube names as defined in cuFFS parset file
    cube_dir: str
        Path to directory with rm cubes
    """

    # Open cube and grab header and data
    with fits.open(f"{cube_dir}/{rm_prefix}p.phi.dirty.fits") as hdus:
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
    _, ra_x, dec_y = wcs.all_world2pix(0, pogs_pos.ra.deg, pogs_pos.dec.deg, 0)

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

    pogs_pos = pogs_obj_loc("POGSII-EG-321", "../slurm/iono/POGS-II_ExGal.fits")
    ra_x, dec_y, phi_z, phi, data_p, wcs = read_rm_cube(
        pogs_pos, "rts_imgr_", "../data"
    )

    plt.style.use("seaborn")

    # Plot RM image of POGS source
    fig = plt.figure(figsize=(6.4, 6))
    ax = fig.add_subplot(1, 1, 1, projection=wcs[:, :, int(phi_z)])
    im = ax.imshow(data_p[:, :, int(phi_z)], origin="lower", cmap="viridis")

    ax.plot(ra_x, dec_y, "#e85d04", marker="o", mfc="none", ms=32, mew=2.0)

    ax.set_xlim(ra_x - 56, ra_x + 56)
    ax.set_ylim(dec_y - 56, dec_y + 56)

    ax.coords.grid(True, color="white", alpha=0.8, ls="dotted")
    ax.coords[0].set_format_unit(u.deg)
    ax.coords[0].set_auto_axislabel(False)
    ax.set_xlabel("Right Ascension [deg]")

    ax.coords[1].set_format_unit(u.deg)
    ax.coords[1].set_auto_axislabel(False)
    ax.set_ylabel("Declination [deg]")

    ax.set_title(f"POGSII-EG-321 RM: {phi[phi_z]} [rad m$^{-2}$]")

    plt.show()

    # Plot RM Spectra of POGS source
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(phi, data_p[dec_y, ra_x, :], label="POGS RM", color="#207561")
    ax.set_xlabel("Faraday Depth [rad/m$^2$]")
    ax.set_ylabel("Polarized Flux Density [Jy/PSF/RMSF]")
    ax.set_title(f"POGSII-EG-321  $\phi_{{max}}={phi[phi_z]}$")

    leg = plt.legend(frameon=True, markerscale=1, handlelength=1)
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    plt.tight_layout()
    #  plt.savefig("rm_spec.png")
    plt.show()

    # Plot RMSF
    file = "../data/rts_imgr_rmsf.txt"
    data = np.loadtxt(file)

    phi = data[:, 0]
    q = data[:, 1]
    u = data[:, 2]
    p = data[:, 3]

    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(phi, p, color="#207561", label=r"$ \vert R \vert $", zorder=2)
    ax.plot(
        phi,
        q,
        alpha=0.7,
        color="#252525",
        linewidth=1.8,
        linestyle=(0, (0.1, 1.2)),
        dash_capstyle="round",
        label=r"$ real(R) $",
        zorder=3,
    )
    ax.plot(
        phi,
        u,
        alpha=0.7,
        color="#350b40",
        linewidth=1.6,
        linestyle="dashed",
        dash_capstyle="round",
        label=r"$ imag(R) $",
        zorder=1,
    )
    ax.set_xlabel("Faraday Depth [rad/m$^2$]")
    ax.set_ylabel("RMSF")
    ax.set_title("POGSII-EG-321 RMSF")

    leg = plt.legend(frameon=True, markerscale=1, handlelength=1)
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    plt.tight_layout()
    #  plt.savefig("rmsf.png")
    plt.show()
