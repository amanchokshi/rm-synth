"""
Determine flux leakage from stokes I to Q, U, V.

GLEAM sources within an image area are weighted by the
FEE beam using mwa_pb.primary_beam and ranked by I flux.
The brightest 1000 are selected and the flux at
corresponding locations in the Q, U, V images are
determined. It is assumed that most extragalactic
sources are unpolarized, or that depolarization over
the large bandwidth of the MFS images depolarizes any
polarised flux. Thus, flux measured in Q, U, V images
can be attributed to instrumental leakage.

A quadratic surface is fit to the fractional flux in
Q, U, V.

Adapted from code by:
Jack Line - https://github.com/JLBLine/srclists
Chris Riseley - correct_leakage.py
"""

from pathlib import Path

import mwa_hyperbeam
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mwa_pb import primary_beam
from scipy.linalg import lstsq
from skyfield.api import wgs84

# Ignore SettingWithCopyWarning
# https://stackoverflow.com/questions/20625582/how-to-deal-with-settingwithcopywarning-in-pandas
pd.options.mode.chained_assignment = None


def read_metafits(metafits):
    """Extract header info from MWA metafits file

    Parameters
    ----------
    metafits : str
        Path to metafits file

    Returns
    -------
    lst : float
        Local Sidereal Time
    obsid : int
        GPS time at start of observation
    freqcent : float
        Centre frequency of observation
    ra_point : float
        Ra phase center of observation in decimal degrees
    dec_point : float
        Dec phase center of observation in decimal degrees
    delays : list
        Beamformer delays
    """

    with fits.open(metafits) as hdus:
        hdu = hdus[0]
        hdr = hdu.header

        lst = float(hdr["LST"])
        obsid = hdr["GPSTIME"]
        freqcent = hdr["FREQCENT"] * 1e6
        ra_point = hdr["RA"]
        dec_point = hdr["DEC"]
        delays = np.fromstring(hdr["DELAYS"], dtype=int, sep=",")

    return lst, obsid, freqcent, ra_point, dec_point, delays


def makeUnpolInstrumentalResponse(j1, j2):
    """Convert Jones matracies to unpolarised beam responses.

    Form the visibility matrix in instrumental response from two Jones
    matrices assuming unpolarised sources (hence the brightness matrix is
    the identity matrix)

    Designed to work with jones matracies output from mwa_hyperbeam, which
    have shape of first axis being equal to number of az/za being evaluated
    with the last axis having 4 values. [[az][za], 4]

    Modified from: https://github.com/MWATelescope/mwa_pb - beam_tools.py


    Parameters
    ----------
    j1 : numpy.array
        Jones matrix output from mwa_pb. [[az], 4]
    j2 : numpy.array
        Jones matrix output from mwa_pb. [[az], 4]

    Returns
    -------
    numpy.array
       Shape [[az][za], [xx, xy, yx, yy]]
    """
    result = np.empty_like(j1)

    result[:, 0] = j1[:, 0] * j2[:, 0].conjugate() + j1[:, 1] * j2[:, 1].conjugate()
    result[:, 3] = j1[:, 2] * j2[:, 2].conjugate() + j1[:, 3] * j2[:, 3].conjugate()
    result[:, 1] = j1[:, 0] * j2[:, 2].conjugate() + j1[:, 1] * j2[:, 3].conjugate()
    result[:, 2] = j1[:, 2] * j2[:, 0].conjugate() + j1[:, 3] * j2[:, 1].conjugate()
    return result


def band_avg_flux(freqcent, gleam_df):
    """Find flux at freqcent with spectral index and nearest gleam band.

    Parameters
    ----------
    freqcent : float/str
        Frequency centre from metafits
    gleam_df : pandas dataframe
        GLEAM catalog as a pandas dataframe

    Returns
    -------
    numpy.array
        Flux at freqcent
    """

    # 200 - 230 MHz
    if float(freqcent) > 200e6:
        gleam_bands = np.array([204, 212, 220, 227])

    # 163 - 200 MHz
    else:
        gleam_bands = np.array([166, 174, 181, 189, 197])

    # Create array of fluxes in gleam bands
    flux_array = []
    for b in gleam_bands:
        flux_array.append(gleam_df[f"peak_flux_{b}"].to_numpy())

    # Average peak flux in each of the gleam bands
    sfluxes = np.sum(np.asarray(flux_array), axis=0) / gleam_bands.shape[0]

    return sfluxes


def eq2horz(HA, Dec, lat):
    """Convert equatorial coordinates to horizontal

    [Az,Alt]=eq2horz(HA,Dec,lat)

    The sign convention for azimuth is north zero, east +pi/2 from
    erfa eraHd2ae https://github.com/liberfa/erfa/blob/master/src/hd2ae.c
    azimuth here is defined with N=0

    Parameters
    ----------
    HA : List or np.array()
        Array of hour angles in decimal degrees
    Dec : List or np.array()
        Array of declinations in decimal degrees
    lat : float
        Latitude of observer in decimal degrees

    Returns
    -------
    Array
        [Az, Alt] : Azimuth and Altitude arrays

    """

    HA = np.asarray(HA)
    Dec = np.asarray(Dec)

    sh = np.sin(HA * np.pi / 180)
    ch = np.cos(HA * np.pi / 180)
    sd = np.sin(Dec * np.pi / 180)
    cd = np.cos(Dec * np.pi / 180)
    sl = np.sin(lat * np.pi / 180)
    cl = np.cos(lat * np.pi / 180)

    # (Az,El) as (x,y,z)
    x = -ch * cd * sl + sd * cl
    y = -sh * cd
    z = ch * cd * cl + sd * sl

    # to spherical
    r = np.sqrt(x * x + y * y)
    a = np.arctan2(y, x)
    a[np.where(r == 0)] = 0
    a[np.where(a < 0)] += np.pi * 2
    el = np.arctan2(z, r)

    # Convert back to degrees
    return [a * 180 / np.pi, el * 180 / np.pi]


def get_beam_weights(
    ras=None,
    decs=None,
    LST=None,
    mwa_lat=None,
    freqcent=None,
    delays=None,
    mwa_pb=False,
):
    """Get normalized beam weights for an mwa observation

    Takes ra and dec coords, and works out the overall
    normalized beam power at that location using the
    2016 spherical harmonic beam code from mwa_pb

    Parameters
    ----------
    ras : array
        Ras at which beam is to be evaluated
    decs : array
        Decs at which beam is to be evaluated
    lst : float
        Local Sidereal Time
    mwa_lat : float
        MWA latitude in decimal degrees
    freqcent : float
        Centre frequency of observation
    delays : list
        Beamformer delays

    Returns
    -------
    Array:
        Beam weights at positions of ras, decs
    """

    # For each gleam source, work out it's position, convolve with the beam and sum for the source

    # HA=LST-RA in def of ephem_utils.py
    has = LST - np.array(ras)

    # Convert ras, decs to az, alt
    Az, Alt = eq2horz(has, np.array(decs), mwa_lat)

    # Convert to zenith angle, azmuth in rad
    za = (90 - Alt) * np.pi / 180
    az = Az * np.pi / 180

    if mwa_pb:

        XX, YY = primary_beam.MWA_Tile_full_EE(
            [za],
            [az],
            freq=freqcent,
            delays=delays,
            zenithnorm=True,
            power=True,
            interp=False,
        )
        # Old way of combining XX and YY - end up with beam values greater than 1, not good!
        # beam_weights = sqrt(XX[0]**2+YY[0]**2)
        beam_weights = np.asarray((XX[0] + YY[0]) / 2.0)

        return beam_weights

    else:

        # We can make a new beam object with a path to the HDF5 file specified by MWA_BEAM_FILE.
        beam = mwa_hyperbeam.FEEBeam()

        amps = [1.0] * 16
        norm_to_zenith = True

        # beam.calc_jones is also available, but that makes a single Jones matrix at a
        # time, so one would need to iterate over az and za. calc_jones_array is done in
        # parallel with Rust (so it's fast).
        jones = beam.calc_jones_array(az, za, freqcent, delays, amps, norm_to_zenith)

        unpol_beam = makeUnpolInstrumentalResponse(jones, jones)

        XX = unpol_beam[:, 0]
        YY = unpol_beam[:, 3]

        beam_weights = np.asarray((XX + YY) / 2.0)

        return beam_weights


def gleam_by_beam(
    gleam_cat=None,
    mfs_fits=None,
    lst=None,
    mwa_lat=None,
    freqcent=None,
    delays=None,
    smin=None,
    beam_thresh=None,
    plot=False,
):
    """Find bright GLEAM sources in sensitive part of the beam.

    Parameters
    ----------
    gleam_cat : str
        Path to GLEAM fits catalog
    mfs_fits : str
        Path to mfs stokes I image
    lst : float
        Local Sidereal Time
    mwa_lat : float
        MWA latitude in decimal degrees
    freqcent : float
        Centre frequency of observation
    delays : list
        Beamformer delays
    smin : float
        Minimum source flux to be considered in Jy
    beam_thresh : float
        Minimum zenith normalized beam power
    plot : boolean
        If true, plot something really cool

    Returns
    -------
    pandas.dataframe :
        gleam_beam - pandas dataframe with gleam sources in sensitive area of beam
        and brighter than minimum flux. Also includes beam_weights and ra_pix, dec_pix
        pixel coordinates of peak flux locations of data array.
    """

    # Read the gleam catalog and return astropy Table object
    dat = Table.read(gleam_cat)

    # convert to a pandas data frame
    gleam_df = dat.to_pandas()

    # Read MFS image and extract range of ras and decs
    with fits.open(mfs_fits) as hdus:
        hdr = hdus[0].header
        data = hdus[0].data

        # World Coordinate System
        wcs = WCS(hdr)

    # wcs.all_pix2world(⍺, δ, freq, stokes, 0)
    ras, _, _, _ = wcs[:2].all_pix2world(np.arange(data.shape[3]), 0, 0, 0, 0)
    _, decs, _, _ = wcs[:2].all_pix2world(0, np.arange(data.shape[2]), 0, 0, 0)

    # Rescale ras from -180, 180 to 0, 360
    ras = np.mod(ras, 360)

    # Select gleam sources within image
    gleam_df = gleam_df[
        (gleam_df.RAJ2000 <= np.amax(ras))
        & (gleam_df.RAJ2000 >= np.amin(ras))
        & (gleam_df.DEJ2000 <= np.amax(decs))
        & (gleam_df.DEJ2000 >= np.amin(decs))
    ]

    # Extract columns of interest
    gleam_in_img = gleam_df[["Name", "RAJ2000", "DEJ2000"]]

    # Add column with band averaged flux
    sfluxes = band_avg_flux(freqcent, gleam_df)
    gleam_in_img["sfluxes"] = sfluxes

    gleam_in_img = gleam_in_img[(gleam_in_img.sfluxes >= smin)]

    # Determine beam weights for gleam_in_img sources
    beam_weights = get_beam_weights(
        ras=gleam_in_img["RAJ2000"].to_numpy(),
        decs=gleam_in_img["DEJ2000"].to_numpy(),
        LST=lst,
        mwa_lat=mwa_lat,
        freqcent=freqcent,
        delays=delays,
    )

    # Append beam weights to df as new column
    gleam_in_img["beam_weights"] = beam_weights

    # Bleam sources in the field with beam weights more than threshold
    gleam_beam = gleam_in_img[(gleam_in_img.beam_weights >= beam_thresh)]

    # Convert gleam ra, dec to pix coordinates
    ra_pix = []
    dec_pix = []
    for i in range(gleam_beam["RAJ2000"].to_numpy().shape[0]):
        coord = wcs.wcs_world2pix(
            gleam_beam["RAJ2000"].to_numpy()[i],
            gleam_beam["DEJ2000"].to_numpy()[i],
            0,
            0,
            0,
        )
        ra_pix.append(coord[0])
        dec_pix.append(coord[1])

    # Round to integer coordinate values
    ra_pix_int = np.rint(np.asarray(ra_pix)).astype(int)
    dec_pix_int = np.rint(np.asarray(dec_pix)).astype(int)

    # Now we check the square of 9 pixels centered around each source
    # ra_pix_int, dec_pix_int to find the index of the peak pixel
    # This is usually the centre pixel, but not always

    ra_pix_peak = []
    dec_pix_peak = []

    for p in range(ra_pix_int.shape[0]):

        # Extract 3x3 pixels centred around coords
        data_2d = data[0, 0, :, :]
        sarray = data_2d[
            dec_pix_int[p] - 1 : dec_pix_int[p] + 2,
            ra_pix_int[p] - 1 : ra_pix_int[p] + 2,
        ]

        # Find indices of peak in subarray
        # https://stackoverflow.com/questions/55284090/how-to-find-maximum-value-in-whole-2d-array-with-indices
        peak_inds = np.unravel_index(sarray.argmax(), sarray.shape)

        # Convert peak_inds of sarray to indices of original data array
        # ra_pix_int[p] is ra index of center of sarray
        # ra_pix_int[p] - 1 take you back to the top left corner of subarray
        ra_pix_peak.append(ra_pix_int[p] - 1 + peak_inds[1])
        dec_pix_peak.append(dec_pix_int[p] - 1 + peak_inds[0])

    # Add pixel coordinates of peak fluxes to data frame
    gleam_beam["ra_pix"] = ra_pix_peak
    gleam_beam["dec_pix"] = dec_pix_peak

    # Plotting stuff
    if plot is True:
        plt.style.use("seaborn")
        fig = plt.figure(figsize=(7, 7))
        ax = fig.add_subplot(111, projection=wcs[0, 0, :, :])
        ax.imshow(data[0, 0, :, :], origin="lower", cmap="Spectral_r")
        ax.scatter(
            ra_pix,
            dec_pix,
            marker="o",
            facecolor="none",
            edgecolor="seagreen",
            linewidth=1.2,
            label="GLEAM RaDec",
        )
        ax.scatter(
            ra_pix_int,
            dec_pix_int,
            marker="o",
            facecolor="none",
            edgecolor="darkorange",
            linewidth=2.1,
            label="GLEAM Pixel",
        )
        ax.scatter(
            ra_pix_peak,
            dec_pix_peak,
            marker="o",
            facecolor="none",
            edgecolor="crimson",
            linewidth=1.2,
            label="Source Peak",
        )
        ax.coords.grid(True, color="white", alpha=0.8, ls="dotted")
        ax.coords[0].set_format_unit(u.deg)
        ax.coords[0].set_auto_axislabel(False)
        ax.set_xlabel("Right Ascension [deg]")

        ax.coords[1].set_format_unit(u.deg)
        ax.coords[1].set_auto_axislabel(False)
        ax.set_ylabel("Declination [deg]")

        ax.set_title(
            f"GLEAM sources - Beam Weight >= {beam_thresh}, Source Flux >= {smin} Jy"
        )

        leg = plt.legend(frameon=True, markerscale=1, handlelength=1)
        leg.get_frame().set_facecolor("white")
        for le in leg.legendHandles:
            le.set_alpha(1)

        plt.show()

    return gleam_beam


def fit_leakage(
    gleam_beam=None,
    mfs_dir=None,
    title=None,
    LST=None,
    mwa_lat=None,
    freqcent=None,
    delays=None,
):
    """Fit leakage surfaces to Q, U, V images."""

    # Ra, Dec pixel positions from stokes I mfs image
    ra_pix = gleam_beam.ra_pix.to_numpy()
    dec_pix = gleam_beam.dec_pix.to_numpy()

    pols = ["I", "Q", "U", "V"]

    for pol in pols:

        msf_fits = f"{mfs_dir}/uvdump-MFS-{pol}-image.fits"

        with fits.open(msf_fits) as hdus:
            data = hdus[0].data

            # Flux at ra, dec pixel positions of mfs images
            flux = data[0, 0, dec_pix, ra_pix]

            # Add observed flux to pandas dataframe
            gleam_beam[f"{pol}_flux"] = flux

    # Determine fractional leakage for Q, U, V

    I_flux = gleam_beam.I_flux.to_numpy()
    Q_flux = gleam_beam.Q_flux.to_numpy()
    U_flux = gleam_beam.U_flux.to_numpy()
    V_flux = gleam_beam.V_flux.to_numpy()

    gleam_beam["Q_leak"] = Q_flux / I_flux
    gleam_beam["U_leak"] = U_flux / I_flux
    gleam_beam["V_leak"] = V_flux / I_flux

    # Get image data from I header and create wcs obj
    msf_i = f"{mfs_dir}/uvdump-MFS-I-image.fits"

    with fits.open(msf_i) as hdus:
        hdr = hdus[0].header

        # World Coordinate System
        wcs = WCS(hdr)

        # Pixel indices along Ra: x, Dec: y
        x = np.arange(hdr["NAXIS1"])
        y = np.arange(hdr["NAXIS2"])

        # Pixel meshgrid
        X, Y = np.meshgrid(x, y)
        ras, decs, _, _ = wcs.wcs_pix2world(X, Y, 0, 0, 0)

        # Rescale ras from -180, 180 to 0, 360
        ras = np.mod(ras, 360)

        # Flatten 2D ras, decs arrays
        ra_f = ras.flatten()
        dec_f = decs.flatten()

    # Fit leakage surfaces for Q, U, V polarizations

    beam_weights = get_beam_weights(
        ras=ra_f, decs=dec_f, LST=lst, mwa_lat=mwa_lat, freqcent=freqcent, delays=delays
    )
    beam_weights = beam_weights.reshape(ras.shape)

    # Dictionary for leakage surface arrays
    leakage_surface = {}

    for i, pol in enumerate(["Q", "U", "V"]):

        # http://inversionlabs.com/2016/03/21/best-fit-surfaces-for-3-dimensional-data.html
        # Create data array to fit quadratic surface to
        # x, y, z in array columns
        data = np.c_[
            gleam_beam.RAJ2000.to_numpy(),
            gleam_beam.DEJ2000.to_numpy(),
            gleam_beam[f"{pol}_leak"].to_numpy(),
        ]

        # Best fit quadratic surface (2nd-order)
        A = np.c_[
            np.ones(data.shape[0]),
            data[:, :2],
            np.prod(data[:, :2], axis=1),
            data[:, :2] ** 2,
        ]
        B = data[:, 2]
        C, _, _, _ = lstsq(A, B)

        # Evaluate the fit on the ra, dec 2D grid
        Z = np.dot(
            np.c_[
                np.ones(ra_f.shape), ra_f, dec_f, ra_f * dec_f, ra_f ** 2, dec_f ** 2
            ],
            C,
        ).reshape(ras.shape)

        leakage_surface[f"{pol}"] = Z

    plt.style.use("seaborn")
    fig, axs = plt.subplots(
        1, 3, figsize=(16, 6), subplot_kw=dict(projection=wcs[0, 0, :, :])
    )
    fig.suptitle(title, fontsize=16)

    for i, pol in enumerate(["Q", "U", "V"]):

        im = axs[i].imshow(leakage_surface[f"{pol}"], origin="lower", cmap="Spectral_r")

        levels = [0.001, 0.01, 0.1, 0.3, 0.6, 0.9]
        CS = axs[i].contour(
            beam_weights.real,
            levels,
            colors="#222222",
            linewidths=0.7,
            linestyles="dotted",
        )
        axs[i].clabel(CS, inline=1, fontsize=7)

        axs[i].scatter(
            gleam_beam.ra_pix.to_numpy(),
            gleam_beam.dec_pix.to_numpy(),
            c=gleam_beam[f"{pol}_leak"].to_numpy(),
            marker="s",
            cmap="Spectral_r",
            s=77,
            edgecolor="black",
            linewidth=0.4,
            vmin=np.amin(leakage_surface[f"{pol}"]),
            vmax=np.amax(leakage_surface[f"{pol}"]),
            zorder=777,
        )

        axs[i].coords.grid(True, color="white", alpha=0.8, ls="dotted")
        axs[i].coords[0].set_format_unit(u.deg)
        axs[i].coords[0].set_auto_axislabel(False)
        axs[i].set_xlabel("Right Ascension [deg]")

        axs[i].coords[1].set_format_unit(u.deg)
        axs[i].coords[1].set_auto_axislabel(False)
        axs[i].set_ylabel("Declination [deg]")

        axs[i].set_title(f"Quadratic Leakage Surface [{pol}/I]", y=1.2)

        divider = make_axes_locatable(axs[i])
        cax = divider.append_axes("top", size="5%", pad=0.1)
        cax.coords[0].set_ticklabel_position("t")
        cax.coords[0].set_axislabel_position("t")
        cax.coords[0].set_axislabel("Fractional Leakage")
        cax.coords[1].set_ticks(
            alpha=0, color="w", size=0, values=[] * u.dimensionless_unscaled
        )
        plt.colorbar(im, cax=cax, orientation="horizontal")

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()

    axs[1].coords[1].set_ticks_visible(False)
    axs[1].coords[1].set_ticklabel_visible(False)
    axs[2].coords[1].set_ticks_visible(False)
    axs[2].coords[1].set_ticklabel_visible(False)

    #  plt.show()
    plt.savefig(f"../data/leakage/{title}.png", bbox_inches="tight", dpi=300)


if __name__ == "__main__":

    gleam_cat = Path("../data/leakage/GLEAM_EGC_v2.fits")

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

        gleam_beam = gleam_by_beam(
            gleam_cat=gleam_cat,
            mfs_fits=f"{mfs_dir}/uvdump-MFS-I-image.fits",
            lst=lst,
            mwa_lat=mwa_lat,
            freqcent=freqcent,
            delays=delays,
            smin=2.0,
            beam_thresh=0.3,
            plot=False,
        )

        gleam_beam = fit_leakage(
            gleam_beam=gleam_beam,
            mfs_dir=mfs_dir,
            title=f"{d}",
            LST=lst,
            mwa_lat=mwa_lat,
            freqcent=freqcent,
            delays=delays,
        )
