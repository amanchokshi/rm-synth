import matplotlib
import numpy as np
from astropy import units as u
from astropy.coordinates import Galactic, SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy_healpix import HEALPix
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from colormaps import spectral


def read_haslam(fits_file):
    """Read the healpix haslam map."""

    with fits.open("haslam408_dsds_Remazeilles2014.fits") as hdus:
        nside = hdus[1].header["NSIDE"]
        order = hdus[1].header["ORDERING"]
        temp = np.ravel((hdus[1].data["TEMPERATURE"][:, :]))

    return temp, nside, order


def eor_field(x0, y0, fov, radec):
    """Plot an EoR field outline.

    x0, y0:     Pixel coords of centre of image
    fov:        Field of view of half power beam
    ra, dec:    Right ascension & Declination list
    """
    return plt.Circle(
        (x0 + 15 * radec[0], y0 + radec[1]),
        fov / 2,
        ec="w",
        lw=2.1,
        ls=":",
        fill=None,
        zorder=4,
    )


if __name__ == "__main__":

    # HASLAM MAP
    # ----------
    # Download Haslam fits file at https://lambda.gsfc.nasa.gov/product/foreground/fg_2014_haslam_408_get.cfm
    temp, nside, order = read_haslam("haslam408_dsds_Remazeilles2014.fits")

    # Read the POGs tables in to return an astropy Table object
    exgal = Table.read("POGS-II_ExGal.fits")
    psr = Table.read("POGS-II_PsrCat.fits")

    df_ex = exgal.to_pandas()
    df_psr = psr.to_pandas()

    # Size of the marker based on absolute value of rm
    s_exgal = np.absolute(df_ex.rm)
    s_psr = np.absolute(df_psr.rm)

    # Sample a 1024x1024 grid in RA/Dec
    ra = np.linspace(-180, 180, 3601) * u.deg
    dec = np.linspace(-90.0, 90.0, 1801) * u.deg
    ra_grid, dec_grid = np.meshgrid(ra, dec)

    # Set up Astropy coordinate objects
    coords = SkyCoord(ra_grid.ravel(), dec_grid.ravel(), frame="icrs")

    # Interpolate values
    hp = HEALPix(nside=nside, order=order, frame=Galactic())
    tmap = hp.interpolate_bilinear_skycoord(coords, temp)
    tmap = tmap.reshape((1801, 3601))

    # EoR Fields - [α, δ]
    # -------------------
    eor0 = [0, -27]
    eor1 = [4, -27]
    eor2 = [10.3, -10]
    fov = 26  # half power beam
    pix_scale = 0.1  # deg/pixel

    # Plotting
    # --------
    nice_fonts = {
        "font.family": "sans-serif",
        "axes.labelsize": 10,
        "axes.titlesize": 10,
        "font.size": 8,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
    }

    plt.rcParams.update(nice_fonts)

    # Custom spectral colormap
    spec, spec_r = spectral()

    # Figure
    fig = plt.figure(figsize=(14, 7))
    ax = fig.add_subplot(1, 1, 1)

    # Haslam map
    haslam = plt.imshow(
        tmap,
        cmap=spec,
        extent=[-180, 180, -90, 90],
        norm=matplotlib.colors.LogNorm(),
        aspect="auto",
        origin="lower",
        zorder=1,
    )

    ex = plt.scatter(
        [i - 360 if i > 180 else i for i in df_ex.ra],
        df_ex.dec,
        s=np.sqrt(s_exgal) * 7,
        #  s=49,
        c=df_ex.rm,
        cmap="RdYlBu_r",
        ec="#222",
        label="ExGal",
        vmin=-200,
        vmax=200,
        alpha=0.9,
        zorder=2,
    )
    plt.scatter(
        [i - 360 if i > 180 else i for i in df_psr.ra],
        df_psr.dec,
        s=np.sqrt(s_psr) * 7,
        #  s=49,
        c=df_psr.rm,
        cmap="RdYlBu_r",
        ec="#222",
        marker="s",
        vmin=min(df_ex.rm),
        vmax=max(df_ex.rm),
        label="PSR",
        alpha=0.9,
        zorder=3,
    )

    # Labels, ticks, etc
    ax.set_xticks(np.linspace(-180, 180, 7))
    ax.set_yticks(np.linspace(-90, 90, 7))
    ax.set_xticklabels(np.arange(-12, 13, 4))
    ax.set_yticklabels(np.arange(-90, 91, 30))
    ax.set_title("Haslam 408$MHz$")
    ax.set_xlabel("Right ascension [$h$]")
    ax.set_ylabel("Declination [$deg$]")

    # Colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2.4%", pad=0.14)
    fig.colorbar(haslam, cax=cax, label="Temperature [$\degree$$K$]")

    #  divider2 = make_axes_locatable(ax)
    cax2 = divider.append_axes("right", size="2.4%", pad=0.60)
    fig.colorbar(ex, cax=cax2, label="RM [$rad\ m^{-2}$]")

    # EoR fields
    e0 = eor_field(0, 0, fov, eor0)
    e1 = eor_field(0, 0, fov, eor1)
    e2 = eor_field(0, 0, fov, eor2)
    ax.add_patch(e0)
    ax.add_patch(e1)
    ax.add_patch(e2)
    ax.axis("scaled")

    plt.tight_layout()
    plt.savefig("haslam_pogs.png")
    #  plt.show()
