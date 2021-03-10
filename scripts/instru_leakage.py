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

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from matplotlib import pyplot as plt
from mwa_pb import primary_beam
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
    ras=None, decs=None, LST=None, mwa_lat=None, freqcent=None, delays=None
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
    print(ra_pix)
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


def fit_leakage(gleam_beam=None, mfs_dir=None):
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

    return gleam_beam


if __name__ == "__main__":

    gleam_cat = Path("../data/leakage/GLEAM_EGC_v2.fits")
    mfs_dir = Path("../data/leakage/1120300352")
    metafits = Path(f"{mfs_dir}/1120300352_metafits_ppds.fits")

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
        plot=True,
    )

    gleam_beam = fit_leakage(gleam_beam=gleam_beam, mfs_dir=mfs_dir)

    plt.scatter(
        gleam_beam.RAJ2000.to_numpy(),
        gleam_beam.DEJ2000.to_numpy(),
        c=gleam_beam.Q_leak.to_numpy(),
        marker="s",
        cmap="viridis",
        s=77,
        edgecolor="black",
        linewidth=0.4,
    )
    plt.tight_layout()
    plt.show()

    # Get ra, dec coords of all pixels
    # https://github.com/astropy/astropy/issues/1587
    #  x = np.arange(NAXIS1)
    #  y = np.arange(NAXIS2)
    #  X, Y = np.meshgrid(x, y)
    #  ra, dec = wcs.wcs_pix2world(X, Y, 0)
