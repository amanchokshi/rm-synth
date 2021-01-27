import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from spectral_cube import SpectralCube


def cube_hdr(fits_dir, prefix, pol, chans):
    """Create FITS header for spectral cube.

    Modify the FITS header of the lowest fine frequency image
    to represent a spectral cube.

    Parameters
    ----------
    fits_dir : str
        Path to data directory with fits images
    prefix : str
        Prefix of fits files
    pol : str
        Stokes polarization - Q/U
    chans : int
        Number of fine channels

    Returns
    -------
    astropy.wcs
        Astropy WCS object for spectral cube
    """

    with fits.open(f"{fits_dir}/{prefix}-0000-{pol}-dirty.fits") as hdul:
        hdr = hdul[0].header
        hdr["NAXIS"] = 3
        hdr["NAXIS3"] = chans
        rm_hdr = ("NAXIS4", "CTYPE4", "CRPIX4", "CRVAL4", "CDELT4", "CUNIT4")
        for i in rm_hdr:
            hdr.pop(i, None)

        return WCS(hdr)


def cube_data(fits_dir, prefix, pol, chans, dim):
    """Create a FITS spectral cube.

    Combine a set of fine channel fits images into a spectral cube.

    Parameters
    ----------
    fits_dir : str
        Path to data directory with fits images
    prefix : str
        Prefix of fits files
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

    # loop through fits files
    for i in range(chans):

        # Path to fits file
        fts = f"{fits_dir}/{prefix}-{i:04}-{pol}-dirty.fits"

        # Read fits file and append contents to data array
        with fits.open(fts) as hdul:
            hdu = hdul[0]
            data = hdu.data[0, 0, :, :]
            cube[i, :, :] = data

    return cube.astype(np.float32)


def create_spec_cube(
    fits_dir=None, prefix=None, pol=None, chans=None, dim=None, out_dir=None
):
    """Creates and saves FITS spectral cube.

    Parameters
    ----------
    fits_dir : str
        Path to data directory with fits images
    prefix : str
        Prefix of fits files
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

    wcs = cube_hdr(fits_dir, prefix, pol, chans)
    data = cube_data(fits_dir, prefix, pol, chans, dim)

    cube = SpectralCube(data=data, wcs=wcs)
    cube.hdu.writeto(f"{out_dir}/cube_{pol}_{chans}.fits")


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
        pol=args.pol,
        chans=args.chans,
        dim=args.dim,
        out_dir=args.out_dir,
    )
