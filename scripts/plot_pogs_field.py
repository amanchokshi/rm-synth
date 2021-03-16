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


def eor_field(fov, radec):
    """Plot an EoR field outline.

    fov:        Field of view of half power beam
    ra, dec:    Right ascension & Declination list
    """
    return plt.Circle(
        (15 * radec[0], radec[1]), fov / 2, ec="w", lw=2.1, ls=":", fill=None, zorder=4,
    )


def plt_field(eor, fov, fname, save=False):

    temp, nside, order = read_haslam("haslam408_dsds_Remazeilles2014.fits")

    # Read the POGs tables in to return an astropy Table object
    exgal = Table.read("POGS-II_ExGal.fits")
    df_ex = exgal.to_pandas()

    # Rescale RA from (0, 360) to (-180, 180)
    df_ex["ra"] = [i - 360 if i > 180 else i for i in df_ex.ra]

    # Crop around EoR Field
    df_ex_cr = df_ex[
        (df_ex.ra < (eor[0] * 15 + (fov / 2)))
        & (df_ex.ra > (eor[0] * 15 - (fov / 2)))
        & (df_ex.dec < (eor[1] + (fov / 2)))
        & (df_ex.dec > (eor[1] - (fov / 2)))
    ]

    # Size of the marker based on absolute value of rm
    s_exgal = np.absolute(df_ex_cr.rm)

    # Sample a grid in RA/Dec
    ra = np.linspace((eor[0] * 15 - (fov / 2)), (eor[0] * 15 + (fov / 2)), 261) * u.deg
    dec = np.linspace((eor[1] - (fov / 2)), (eor[1] + (fov / 2)), 261) * u.deg
    ra_grid, dec_grid = np.meshgrid(ra, dec)

    # Set up Astropy coordinate objects
    coords = SkyCoord(ra_grid.ravel(), dec_grid.ravel(), frame="icrs")

    # Interpolate values
    hp = HEALPix(nside=nside, order=order, frame=Galactic())
    tmap = hp.interpolate_bilinear_skycoord(coords, temp)
    tmap = tmap.reshape((261, 261))

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
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(1, 1, 1)

    # Haslam map
    haslam = plt.imshow(
        tmap,
        cmap=spec,
        extent=[
            (eor[0] * 15 - (fov / 2)),
            (eor[0] * 15 + (fov / 2)),
            (eor[1] - (fov / 2)),
            (eor[1] + (fov / 2)),
        ],
        #  norm=matplotlib.colors.LogNorm(),
        aspect="auto",
        origin="lower",
        zorder=1,
    )

    ras_cr = df_ex_cr.ra.to_numpy()
    decs_cr = df_ex_cr.dec.to_numpy()
    ids_cr = df_ex_cr.catalog_id.to_numpy()

    ex = plt.scatter(
        ras_cr,
        decs_cr,
        s=s_exgal * 7,
        c=df_ex_cr.rm,
        cmap="RdYlBu_r",
        ec="#222",
        label="ExGal",
        zorder=2,
    )

    for i in range(len(ids_cr)):
        plt.annotate(
            ids_cr[i].decode("utf-8"),
            xy=(ras_cr[i], decs_cr[i]),
            xytext=(ras_cr[i] + 0.7, decs_cr[i] - 0.7),
            fontsize=6,
            #  arrowprops=dict(facecolor="#222", arrowstyle="->")
        )

    # Labels, ticks, etc
    ax.set_title(f"{fname} Field $&$ POGS")
    ax.set_xlabel("Right ascension [$deg$]")
    ax.set_ylabel("Declination [$deg$]")

    # Colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2.4%", pad=0.4)
    fig.colorbar(haslam, cax=cax, label="Temperature [$\degree$$K$]")

    cax2 = divider.append_axes("right", size="2.4%", pad=0.6)
    fig.colorbar(ex, cax=cax2, label="RM [$rad\ m^{-2}$]")

    # EoR fields
    e0 = eor_field(fov, eor)
    ax.add_patch(e0)
    ax.axis("scaled")
    plt.tight_layout()

    if save:
        plt.savefig(f"{fname}_pogs.png")
    else:
        plt.show()


if __name__ == "__main__":

    # EoR Fields - [α, δ]
    # -------------------
    eor0 = [0, -27]
    eor1 = [4, -27]
    eor2 = [10.3, -10]
    fov = 26  # half power beam

    plt_field(eor0, fov, "EoR0", save=True)
    plt_field(eor1, fov, "EoR1", save=True)
    plt_field(eor2, fov, "EoR2", save=True)
