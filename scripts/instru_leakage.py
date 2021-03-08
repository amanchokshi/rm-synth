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
        delays = hdr["DELAYS"]

    return lst, obsid, freqcent, ra_point, dec_point, delays


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
    beam_weights = (XX[0] + YY[0]) / 2.0

    return beam_weights


def gleam_by_beam(
    gleam_cat=None, mfs_fits=None, lst=None, mwa_lat=None, freqcent=None, delays=None
):

    # Read the gleam catalog and return astropy Table object
    dat = Table.read(gleam_cat)

    # convert to a pandas data frame
    df = dat.to_pandas()

    # Extract columns of interest
    gleam = df[
        [
            "RAJ2000",
            "DEJ2000",
            "int_flux_166",
            "a_166",
            "b_166",
            "int_flux_174",
            "a_174",
            "b_174",
            "int_flux_181",
            "a_181",
            "b_181",
            "int_flux_189",
            "a_189",
            "b_189",
            "int_flux_197",
            "a_197",
            "b_197",
            "int_flux_204",
            "a_204",
            "b_204",
            "int_flux_212",
            "a_212",
            "b_212",
            "int_flux_220",
            "a_220",
            "b_220",
            "int_flux_227",
            "a_227",
            "b_227",
            "alpha",
        ]
    ]

    # Read MFS image and extract range of ras and decs
    with fits.open(mfs_fits) as hdus:
        hdr = hdus[0].header
        data = hdus[0].data

        # World Coordinate System
        wcs = WCS(hdr)

    # wcs.all_pix2world(⍺, δ, freq, stokes, 0)
    ras, _, _, _ = wcs[:2].all_pix2world(np.arange(data.shape[3]), 0, 0, 0, 0)
    _, decs, _, _ = wcs[:2].all_pix2world(0, np.arange(data.shape[2]), 0, 0, 0)

    # Rescale from -180, 180 to 0, 360
    ras = np.mod(ras, 360)

    # Select gleam sources within image
    gleam_in_img = gleam[
        (gleam.RAJ2000 <= np.amax(ras))
        & (gleam.RAJ2000 >= np.amin(ras))
        & (gleam.DEJ2000 <= np.amax(decs))
        & (gleam.DEJ2000 >= np.amin(decs))
    ]

    beam_weights = get_beam_weights(
        ras=gleam_in_img["RAJ2000"].to_numpy(),
        decs=gleam_in_img["DEJ2000"].to_numpy(),
        LST=lst,
        mwa_lat=mwa_lat,
        freqcent=freqcent,
        delays=delays,
    )

    print(beam_weights)

    # Convert ra, dec to pix coordinates
    ra_pix = []
    dec_pix = []
    for i in range(gleam_in_img["RAJ2000"].to_numpy().shape[0]):
        coord = wcs.wcs_world2pix(
            gleam_in_img["RAJ2000"].to_numpy()[i],
            gleam_in_img["DEJ2000"].to_numpy()[i],
            0,
            0,
            0,
        )
        ra_pix.append(coord[0])
        dec_pix.append(coord[1])

    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection=wcs[0, 0, :, :])
    ax.imshow(data[0, 0, :, :], origin="lower", cmap=plt.cm.viridis)
    ax.scatter(ra_pix, dec_pix, marker="o", facecolor="none", edgecolor="seagreen")
    ax.coords.grid(True, color="white", alpha=0.8, ls="dotted")
    ax.coords[0].set_format_unit(u.deg)
    ax.coords[0].set_auto_axislabel(False)
    ax.set_xlabel("Right Ascension [deg]")

    ax.coords[1].set_format_unit(u.deg)
    ax.coords[1].set_auto_axislabel(False)
    ax.set_ylabel("Declination [deg]")
    plt.show()


if __name__ == "__main__":

    gleam_cat = Path("../data/leakage/GLEAM_EGC_v2.fits")
    metafits = Path("../data/leakage/1120300352_metafits_ppds.fits")

    mfs_i = Path("../data/leakage/wsc_stokes_2048/uvdump-MFS-I-image.fits")

    # MWA coordinates
    mwa_loc = wgs84.latlon(-26.703319, 116.670815, 337.83)
    mwa_lat = mwa_loc.latitude.degrees
    mwa_lon = mwa_loc.longitude.degrees
    mwa_el = mwa_loc.elevation.m

    lst, obsid, freqcent, ra_point, dec_point, delays = read_metafits(metafits)
    gleam_by_beam(
        gleam_cat=gleam_cat,
        mfs_fits=mfs_i,
        lst=lst,
        mwa_lat=mwa_lat,
        freqcent=freqcent,
        delays=delays,
    )
