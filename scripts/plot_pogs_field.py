from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import Galactic, SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy_healpix import HEALPix
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from skyfield.api import wgs84

from instru_leakage import get_beam_weights, read_metafits


def read_haslam(haslam_fits):
    """Read the healpix haslam map."""

    with fits.open(haslam_fits) as hdus:
        nside = hdus[1].header["NSIDE"]
        order = hdus[1].header["ORDERING"]
        temp = np.ravel((hdus[1].data["TEMPERATURE"][:, :]))

    return temp, nside, order


def plt_field(ra_point, dec_point, fov, haslam_fits, pogs_fits):

    temp, nside, order = read_haslam(haslam_fits)

    # Read the POGs tables in to return an astropy Table object
    exgal = Table.read(pogs_fits)
    df_ex = exgal.to_pandas()

    # Rescale RA from (0, 360) to (-180, 180)
    #  df_ex["ra"] = [i - 360 if i > 180 else i for i in df_ex.ra]

    # Crop around EoR Field
    df_ex_cr = df_ex[
        (df_ex.ra < (ra_point + (fov / 2)))
        & (df_ex.ra > (ra_point - (fov / 2)))
        & (df_ex.dec < (dec_point + (fov / 2)))
        & (df_ex.dec > (dec_point - (fov / 2)))
    ]

    # Size of the marker based on absolute value of rm
    s_exgal = np.absolute(df_ex_cr.rm)

    # Sample a grid in RA/Dec
    ra = np.linspace((ra_point - (fov / 2)), (ra_point + (fov / 2)), 261) * u.deg
    dec = np.linspace((dec_point - (fov / 2)), (dec_point + (fov / 2)), 261) * u.deg

    # Setup ra/dec grid
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

    # Figure
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(1, 1, 1)

    # Haslam map
    haslam = plt.imshow(
        tmap,
        cmap="Spectral_r",
        extent=[
            (ra_point - (fov / 2)),
            (ra_point + (fov / 2)),
            (dec_point - (fov / 2)),
            (dec_point + (fov / 2)),
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
    ax.set_xlabel("Right ascension [$deg$]")
    ax.set_ylabel("Declination [$deg$]")

    # Colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2.4%", pad=0.4)
    fig.colorbar(haslam, cax=cax, label="Temperature [$\degree$$K$]")

    cax2 = divider.append_axes("right", size="2.4%", pad=0.6)
    fig.colorbar(ex, cax=cax2, label="RM [$rad\ m^{-2}$]")

    # EoR fields
    #  e0 = eor_field(fov, eor)
    #  ax.add_patch(e0)
    #  ax.axis("scaled")
    plt.tight_layout()

    plt.show()


if __name__ == "__main__":

    fov = 26  # half power beam
    haslam_fits = "../data/leakage/haslam408_dsds_Remazeilles2014.fits"
    pogs_fits = "../data/leakage/POGS-II_ExGal.fits"
    dirs = ["fee_1120300352", "fee_1120300232", "ana_1120300352", "ana_1120300232"]

    for d in dirs:

        print(f" ** INFO: Crunching data in - {d}")

        _, obsid = d.split("_")

        mfs_dir = Path(f"../data/leakage/{d}")
        metafits = Path(f"{mfs_dir}/{obsid}_metafits_ppds.fits")

        # MWA coordinates
        mwa_loc = wgs84.latlon(-26.703319, 116.670815, 337.83)
        mwa_lat = mwa_loc.latitude.degrees
        mwa_lon = mwa_loc.longitude.degrees
        mwa_el = mwa_loc.elevation.m

        lst, obsid, freqcent, ra_point, dec_point, delays = read_metafits(metafits)

        plt_field(ra_point, dec_point, fov, haslam_fits, pogs_fits)

        break
