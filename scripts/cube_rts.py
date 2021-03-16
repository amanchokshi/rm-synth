"""
Convert RTS fine channel images in spectral cube.
"""

import re
from pathlib import Path

import numpy as np
from astropy.io import fits


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


def cube_hdr(fits_dir, pol, freqs):
    """Create FITS header for spectral cube.

    Modify the FITS header of the lowest fine frequency image
    to represent a spectral cube. Make the first axis frequency.

    Parameters
    ----------
    fits_dir : str
        Path to data directory with fits images
    pol : str
        Stokes polarization - Q/U
    freqs : list
        List of frequencies

    Returns
    -------
    astropy.wcs
        Astropy FITS header object for spectral cube
    """

    fits_file = list(Path(fits_dir).glob(f"*_{freqs[0]}MHz_*{pol}.fits"))[0]
    # Copy of the header
    with fits.open(fits_file) as hdul:
        h2 = hdul[0].header

    # Modify the header to have frequency as first axis
    with fits.open(fits_file) as hdul:
        hdr = hdul[0].header
        hdr["NAXIS"] = 3
        hdr["NAXIS3"] = " "

        hdr["NAXIS1"] = " "
        hdr["NAXIS2"] = int(h2["NAXIS1"])
        hdr["NAXIS3"] = int(h2["NAXIS2"])

        hdr["CRPIX1"] = " "
        hdr["CRPIX2"] = h2["CRPIX1"]
        hdr["CRPIX3"] = h2["CRPIX2"]

        hdr["CDELT1"] = " "
        hdr["CDELT2"] = h2["CDELT1"]
        hdr["CDELT3"] = h2["CDELT2"]

        #  hdr["CUNIT1"] = " "
        #  hdr["CUNIT2"] = h2["CUNIT1"]
        #  hdr["CUNIT3"] = h2["CUNIT2"]

        hdr["CTYPE1"] = " "
        hdr["CTYPE2"] = h2["CTYPE1"]
        hdr["CTYPE3"] = h2["CTYPE2"]

        hdr["CRVAL1"] = " "
        hdr["CRVAL2"] = h2["CRVAL1"]
        hdr["CRVAL3"] = h2["CRVAL2"]

        # Remove 4th dimension stuff - Stokes info
        rm_hdr = (
            "FREQ",
            "PV2_1",
            "PV2_2",
            "CD1_1",
            "CD1_2",
            "CD2_1",
            "CD2_2",
            "BZERO",
            "SIZEX",
            "SIZEY",
        )
        for i in rm_hdr:
            hdr.pop(i, None)

        return hdr


def cube_data(fits_dir, pol, freqs, hdr):
    """Create a FITS spectral cube.

    Combine a set of fine channel fits images into a spectral cube.

    Parameters
    ----------
    fits_dir : str
        Path to data directory with fits images
    prefix : str
        Prefix of fits files
    suffix : str
        Suffix of fits files
    pol : str
        Stokes polarization - Q/U
    freqs : list
        List of frequencies
    hdr : astropy
        Astropy header object retured by cube_hdr()

    Returns
    -------
    numpy.array
        Numpy spectral cube with shape [chans, dim, dim]
    """

    dim = int(hdr["NAXIS2"])
    chans = len(freqs)

    cube = np.zeros((chans, dim, dim))

    # loop through fits files
    for i, freq in enumerate(freqs):

        # Path to fits file
        fts = list(Path(fits_dir).glob(f"*_{freq}MHz_*{pol}.fits"))[0]

        # Read fits file and append contents to data array
        with fits.open(fts) as hdul:
            hdu = hdul[0]
            data = hdu.data
            cube[i] = np.copy(data)

    # Intial shape [freq, Dec, Ra]
    # Change to [Dec, Ra, freq]
    # Making NAXIS1 = Freq
    cube = np.moveaxis(cube, 0, -1)

    return cube.astype(np.float32)


def create_spec_cube(
    fits_dir=None, pol=None, out_dir=None,
):
    """Creates and saves FITS spectral cube.

    Parameters
    ----------
    fits_dir : str
        Path to data directory with fits images
    pol : str
        Stokes polarization - Q/U
    freqs : list
        List of frequencies
    hdr : astropy
        Astropy header object retured by cube_hdr()
    out_dir : str
        Path to output directory

    Returns
    -------
        Save FITS spectral cube to out_dir
    """

    freqs = get_freqs(fits_dir)

    # Make out_dir if doesn't exist
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    hdr = cube_hdr(fits_dir, pol, freqs)
    data = cube_data(fits_dir, pol, freqs, hdr)

    # Write FITS cube
    hdul = fits.PrimaryHDU(data=data, header=hdr)
    hdul.writeto(f"{out_dir}/cube_{pol}.fits")

    # Write frequency list if it doesn't exist
    freq_file = Path(f"{out_dir}/frequency.txt")

    if not freq_file.is_file():
        with open(freq_file, "w") as f:
            for line in freqs:
                f.write(f"{line*1e6:.1f}\n")


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Create a FITS spectral cube from a set of fine channel images",
    )

    parser.add_argument(
        "--fits_dir",
        metavar="\b",
        type=str,
        required=True,
        help="Path to directory with fine channel stokes fits images",
    )

    parser.add_argument(
        "--out_dir",
        metavar="\b",
        type=str,
        required=True,
        help="Path to output directory",
    )

    args = parser.parse_args()

    fits_dir = Path(args.fits_dir)

    # Create Q cube
    create_spec_cube(
        fits_dir=fits_dir, pol="Q", out_dir=args.out_dir,
    )

    # Create U cube
    create_spec_cube(
        fits_dir=fits_dir, pol="U", out_dir=args.out_dir,
    )
