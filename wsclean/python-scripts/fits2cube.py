import numpy as np
from astropy.io import fits


def cube_hdr(fits_dir, prefix, suffix, pol, chans):
    """Create FITS header for spectral cube.

    Modify the FITS header of the lowest fine frequency image
    to represent a spectral cube. Make the first axis frequency.

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
    chans : int
        Number of fine channels

    Returns
    -------
    astropy.wcs
        Astropy FITS header object for spectral cube
    """

    # Copy of the header
    with fits.open(f"{fits_dir}/{prefix}-0000-{pol}-{suffix}.fits") as hdul:
        h2 = hdul[0].header

    # Modify the header to have frequency as first axis
    with fits.open(f"{fits_dir}/{prefix}-0000-{pol}-{suffix}.fits") as hdul:
        hdr = hdul[0].header
        hdr["NAXIS"] = 3
        hdr["NAXIS3"] = chans

        hdr["NAXIS1"] = int(h2["NAXIS3"])
        hdr["NAXIS2"] = int(h2["NAXIS1"])
        hdr["NAXIS3"] = int(h2["NAXIS2"])

        hdr["CRPIX1"] = h2["CRPIX3"]
        hdr["CRPIX2"] = h2["CRPIX1"]
        hdr["CRPIX3"] = h2["CRPIX2"]

        hdr["CDELT1"] = h2["CDELT3"]
        hdr["CDELT2"] = h2["CDELT1"]
        hdr["CDELT3"] = h2["CDELT2"]

        hdr["CUNIT1"] = h2["CUNIT3"]
        hdr["CUNIT2"] = h2["CUNIT1"]
        hdr["CUNIT3"] = h2["CUNIT2"]

        hdr["CTYPE1"] = h2["CTYPE3"]
        hdr["CTYPE2"] = h2["CTYPE1"]
        hdr["CTYPE3"] = h2["CTYPE2"]

        hdr["CRVAL1"] = h2["CRVAL3"]
        hdr["CRVAL2"] = h2["CRVAL1"]
        hdr["CRVAL3"] = h2["CRVAL2"]

        # Remove 4th dimension stuff - Stokes info
        rm_hdr = ("NAXIS4", "CTYPE4", "CRPIX4", "CRVAL4", "CDELT4", "CUNIT4")
        for i in rm_hdr:
            hdr.pop(i, None)

        return hdr


def cube_data(fits_dir, prefix, suffix, pol, chans, dim):
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
    chans : int
        Number of fine channels
    dim : int
        Dimensions of image data

    Returns
    -------
    numpy.array
        Numpy spectral cube with shape [chans, dim, dim]
    """

    cube = np.zeros((chans, dim, dim))

    freqs = []

    # loop through fits files
    for i in range(chans):

        # Path to fits file
        fts = f"{fits_dir}/{prefix}-{i:04}-{pol}-{suffix}.fits"

        # Read fits file and append contents to data array
        with fits.open(fts) as hdul:
            hdu = hdul[0]
            freq = hdu.header["CRVAL3"]
            data = hdu.data[0, 0, :, :]
            cube[i, :, :] = data
            freqs.append(freq)

    # Intial shape [freq, Dec, Ra]
    # Change to [Dec, Ra, freq]
    # Making NAXIS1 = Freq
    cube = np.moveaxis(cube, 0, -1)

    return cube.astype(np.float32), freqs


def create_spec_cube(
    fits_dir=None, prefix=None, suffix=None, pol=None, chans=None, dim=None, out_dir=None
):
    """Creates and saves FITS spectral cube.

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
    chans : int
        Number of fine channels
    dim : int
        Dimensions of image data
    out_dir : str
        Path to output directory

    Returns
    -------
        Save FITS spectral cube to out_dir
    """

    hdr = cube_hdr(fits_dir, prefix, suffix, pol, chans)
    data, freqs = cube_data(fits_dir, prefix, suffix, pol, chans, dim)

    # Write FITS cube
    hdul = fits.PrimaryHDU(data=data, header=hdr)
    hdul.writeto(f"{out_dir}/cube_{pol}_{suffix}.fits")

    # Write frequency list
    with open(f"{out_dir}/frequency.txt", "w") as f:
        for fr in freqs:
            f.write(f"{fr:.0f}\n")


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
        help="Path to directory with fine channel fits images",
    )

    parser.add_argument(
        "--prefix",
        metavar="\b",
        type=str,
        default="uvdump",
        help="Prefix to fine channel file names. Default: uvdump",
    )

    parser.add_argument(
        "--suffix",
        metavar="\b",
        type=str,
        default="dirty",
        help="Suffix to fine channel file names - dirty/image. Default: dirty",
    )

    parser.add_argument(
        "--pol",
        metavar="\b",
        type=str,
        required=True,
        help="Stokes polarization. Either Q/U",
    )

    parser.add_argument(
        "--chans",
        metavar="\b",
        type=int,
        required=True,
        help="Number of fine channel images to merge",
    )

    parser.add_argument(
        "--dim",
        metavar="\b",
        type=int,
        required=True,
        help="Pixel dimensions of image data",
    )

    parser.add_argument(
        "--out_dir",
        metavar="\b",
        type=str,
        required=True,
        help="Path to output directory",
    )

    args = parser.parse_args()

    create_spec_cube(
        fits_dir=args.fits_dir,
        prefix=args.prefix,
        suffix=args.suffix,
        pol=args.pol,
        chans=args.chans,
        dim=args.dim,
        out_dir=args.out_dir,
    )
