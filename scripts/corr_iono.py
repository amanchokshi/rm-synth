"""
Determine ionospheric RM and correct it.

Ionospheric RM determined using GPS & EMM data at Ra, Dec
of POGS source, at obsid gps time. Rotation is frequency
dependant and detoration is applied to individual fine
channel stokes images.

# TODO: Figure origin of the factor of -2 in correct_iono

Inspired by code by Chris Riseley
"""

import re
from shutil import copyfile as cp

import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from RMextract import getRM as gt


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
    list
        Ra, Dec of POGS object in radians
    """

    exgal = Table.read(pogs_path)
    df_eg = exgal.to_pandas()
    df_eg["catalog_id"] = df_eg["catalog_id"].str.decode("utf-8")
    df_eg.set_index("catalog_id", inplace=True)

    pogs_source = df_eg.loc[f"{pogs_obj}"]
    pogs_pos = SkyCoord(pogs_source.ra, pogs_source.dec, unit="deg")
    pointing = [pogs_pos.ra.rad, pogs_pos.dec.rad]

    return pointing


def iono_rm_obs(obsid, del_t, pointing):
    """Use RMExtract to compute ionospheric rotation with
    the Enhanced Magnetic Model (EMM).

    Parameters
    ----------
    obsid : srt
        GPS time at beginning of a 2 min MWA observation
    del_t: int
        Time interval in seconds before and after obsid to compute RM
    pointing: list
        Ra, Dec of source in radians

    Returns
    -------
    float
        ionRM - Median ionospheric RM between obsid-del_t and obsid+del_t
    """

    # Hardcoded location of the MWA center
    mwa = EarthLocation.from_geodetic(116.670815, -26.703319, 337.83).itrs

    # MWA center, ITRF xyz in meters
    mwa_cen = [mwa.x.value, mwa.y.value, mwa.z.value]

    # getRM still wants MJD time in seconds
    # from within a minute on either side of the obs
    obs = Time(f"{obsid}", format="gps", scale="utc")
    starttime = obs.mjd * 24 * 3600 - del_t
    endtime = obs.mjd * 24 * 3600 + 120 + del_t

    # RM Extract magic
    RMdict = gt.getRM(
        ionexPath="./ionex-data/",
        radec=pointing,
        timestep=30,
        timerange=[starttime, endtime],
        stat_positions=[mwa_cen],
        useEMM=True,
    )

    ionRM = np.median(np.array(RMdict["RM"]["st1"]))

    return ionRM


def get_freqs(fits_dir):
    """Get list of frequencies in MHz from files in fits_dir

    Parameters
    ----------
    fits_dir : pathlib.Path
        Path to directory containing fits files

    Returns
    -------
    list
        Frequencies in MHz, ascending order
    """

    file_list = [rf"{f.stem}" for f in fits_dir.glob("*I.fits") if f.is_file]

    freqs = []
    for f in file_list:
        match = re.search(r"\d{,3}\.\d{,3}MHz", f).group(0)
        freq = re.search(r"\d{,3}\.\d{,3}", match).group(0)
        freqs.append(float(freq))

    return sorted(freqs)


def correct_iono(freq, ionRM, fits_dir, out_dir):

    # Find fits files
    ifn = list(fits_dir.glob(f"*{freq}MHz*I.fits"))[0]
    qfn = list(fits_dir.glob(f"*{freq}MHz*Q.fits"))[0]
    ufn = list(fits_dir.glob(f"*{freq}MHz*U.fits"))[0]
    vfn = list(fits_dir.glob(f"*{freq}MHz*V.fits"))[0]

    # Speed of light
    c = 299792458.0

    # Lambda squared
    lam2 = (c / (freq * 1.0e6)) ** 2.0

    # Angle to derotate by
    # Why -2???
    angle = -2 * ionRM * lam2

    hduq = fits.open(qfn)
    qheader = hduq[0].header

    hduu = fits.open(ufn)
    uheader = hduu[0].header

    nq = np.cos(angle) * hduq[0].data - np.sin(angle) * hduu[0].data
    nu = np.sin(angle) * hduq[0].data + np.cos(angle) * hduu[0].data

    fits.writeto(f"{out_dir}/{Path(qfn).name}", nq, qheader, overwrite=True)
    fits.writeto(f"{out_dir}/{Path(ufn).name}", nu, uheader, overwrite=True)

    # Nothing to correct for i and v so just make a copy for consistency
    cp(ifn, f"{out_dir}/{Path(ifn).name}")
    cp(vfn, f"{out_dir}/{Path(vfn).name}")


if __name__ == "__main__":

    import argparse
    from pathlib import Path

    #####################################################################
    #                                                                   #
    #                            Get arguments                          #
    #                                                                   #
    #####################################################################

    parser = argparse.ArgumentParser(description="Derotate ionospheric RM")

    parser.add_argument(
        "--obsid", metavar="\b", default=1061316296, help="Obsid. Default: 1061316296",
    )

    parser.add_argument(
        "--pogs_path",
        metavar="\b",
        default="./POGS-II_ExGal.fits",
        type=str,
        help="Path to POGS EXGAL fits catalogue. Default: ./POGS-II_ExGal.fits",
    )

    parser.add_argument(
        "--pogs_obj",
        metavar="\b",
        default="POGSII-EG-024",
        help="POGS identifier. Default: POGSII-EG-024",
    )

    parser.add_argument(
        "--del_t",
        metavar="\b",
        default=60,
        type=int,
        help="Time on either side of obs to average iono rm over. Default: 60s",
    )

    parser.add_argument(
        "--fits_dir",
        metavar="\b",
        required=True,
        help="Path to dir with stokes fits images",
    )

    parser.add_argument(
        "--out_dir",
        metavar="\b",
        required=True,
        help="Output directory to saved de-rotated stokes fits images to",
    )

    args = parser.parse_args()

    obsid = args.obsid
    pogs_obj = args.pogs_obj
    del_t = args.del_t
    pogs_path = Path(args.pogs_path)
    fits_dir = Path(args.fits_dir)
    out_dir = Path(args.out_dir)

    #####################################################################
    #                                                                   #
    #                 Apply Frequency dependant ionospheric             #
    #                     RM de-rotation to all files                   #
    #                                                                   #
    #####################################################################

    # Make output dir if doesn't exist
    out_dir.mkdir(parents=True, exist_ok=True)

    # Ra, Dec of pogs_obj
    pointing = pogs_obj_loc(pogs_obj, pogs_path)

    # Black magic
    ionRM = iono_rm_obs(obsid, del_t, pointing)

    # All frequencies from files in fits_dir
    freqs = get_freqs(fits_dir)

    for freq in freqs:
        correct_iono(freq, ionRM, fits_dir, out_dir)
