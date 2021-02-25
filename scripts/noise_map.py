# This program is intended to generate a noise map for an RM cube, by fitting a Rayleigh distribution
# to the wings of the RM spectra where there is little/no signal expected.

# Right now, it does so on a pixel-by-pixel basis. It's not clear yet if it would be better to
# smooth things over several pixels, maybe the beam size. Significant variations in the fit
# over single pixels is unrealistic...

from pathlib import Path

import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit


def rayleigh(x, sigma):
    return x / sigma ** 2 * np.exp(-(x ** 2) / 2.0 / sigma ** 2)


def gaussian(x, sigma):
    return 1.0 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-(x ** 2) / (2 * sigma ** 2))


# At one point I tested fitting Gaussians to Stokes Q/U instead of a Rayleigh to I, because Jamie Farnes
# and I had an arugment over which was better. Both were equivalent, as I expected. (Take that, Jamie!)


def noise_map():
    """Add Docsting once script is done"""

    p_cube = Path("../data/noise_321/fee_imgr_p.phi.dirty.fits")
    print(f"Processing {p_cube}")

    with fits.open(p_cube, memmap=True, mode="denywrite") as hdus:
        hdu = hdus[0]
        hdr = hdu.header
        data = hdu.data

        phi_range = np.arange(
            hdr["CRVAL1"] - (hdr["CRPIX1"] - 1) * hdr["CDELT1"],
            hdr["CRVAL1"] - (hdr["CRPIX1"] - 1 - hdr["NAXIS1"]) * hdr["CDELT1"],
            hdr["CDELT1"],
        )

    # Define exclusion region (where real polarization can occur)
    phi_min = -20
    phi_max = 20

    # Array of phi values to exclude
    phi_exc = np.arange(phi_min, phi_max + hdr["CDELT1"], hdr["CDELT1"])

    # Mask array with False at indicies of phi_exc
    phi_mask = np.array([False if p in phi_exc else True for p in phi_range])

    # Data array from rm cube fits file has order: [δ, ⍺, 𝜙]
    # With N1: 𝜙, N2: ⍺, N3: δ in the header
    # We're now going to make a noise array with data order:[⍺, δ]
    # This way, for the new header, N2 can remain ⍺, while N1 is changed to δ and N3 deleted

    noise = np.zeros([data.shape[1], data.shape[0]])

    # Loop over RA
    for i in range(data.shape[1]):

        # Loop over Dec
        for j in range(data.shape[0]):

            print(f"Ra: {i}/{data.shape[1]}, Dec ind:{j}/{data.shape[0]}")

            # Extract RM Spectra for RA, Dec
            rm_spec = data[j, i, :]

            # Mask |𝜙| < 20 where leakage dominates
            rm_masked = rm_spec[phi_mask]

            # For pixels in the primary beam nulls
            if np.sum(rm_masked) == 0:
                noise[i, j] = 0
                continue

            hist, edges = np.histogram(
                rm_masked, range=(np.amin(rm_masked), np.amax(rm_masked)), bins=30
            )

            # Bin centers
            bin_c = edges[:-1] + 0.5 * (edges[1] - edges[0])
            density = hist / (rm_masked.size * (edges[1] - edges[0]))
            sigma = np.sqrt(np.fmax(hist, 1) * (1 - hist / rm_masked.size)) / (
                rm_masked.size * (edges[1] - edges[0])
            )
            # A note on this last line: I calculate the error on the bin values using Poisson statistics,
            # but I add a 'robustness' correction (probably not rigourously correct) where if a bin has zero values
            # I assign it an error equivalent to if it had 1 entry.
            # This prevents bins from having an error of zero, which screws up the fit.

            #  plt.errorbar(bin_c, density, yerr=sigma, fmt='bo', linewidth=3)

            try:
                popt, pcov = curve_fit(rayleigh, bin_c, density, sigma=(1 / sigma))
            except Exception:
                popt = [0]
            # If the fit fails, I give it an error of zero. I haven't seen this happen in any of my maps,
            # so I think it is working reliably in that regard.
            # print "  Noise: ", popt[0],np.sqrt(pcov[0])
            #  print(popt)
            noise[i, j] = np.abs(popt[0])

            #  plt.style.use("seaborn")
            #  plt.plot(bin_c, density, "o-", label="Histogram")
            #  plt.plot(bin_c, rayleigh(bin_c, *popt), label="Rayleigh")
            #  plt.xlabel("RM bins [rad/m^2]")
            #  plt.title("RM Spectra Noise Histogram")
            #  plt.ylabel("Count")

            #  leg = plt.legend(frameon=True, markerscale=1, handlelength=1)
            #  leg.get_frame().set_facecolor("white")
            #  for le in leg.legendHandles:
                #  le.set_alpha(1)

            #  plt.tight_layout()
            #  plt.show()

            #  break
        #  break

    # Create FITS header for noise map, save map to file.
    print(" ** INFO: Writing output noise image")
    new_header = hdr
    new_header["NAXIS"] = 2
    new_header.set("CRVAL1", hdr["CRVAL3"])
    new_header.set("CDELT1", hdr["CDELT3"])
    new_header.set("CRPIX1", hdr["CRPIX3"])
    new_header.set("CTYPE1", hdr["CTYPE3"])
    new_header.set("NAXIS1", hdr["NAXIS3"])
    del new_header["CTYPE3"]
    del new_header["CRVAL3"]
    del new_header["CDELT3"]
    del new_header["CRPIX3"]
    del new_header["NAXIS3"]

    hdul = fits.PrimaryHDU(data=noise, header=new_header)
    hdul.writeto("../data/noise_321/fee_cube_noise.fits")

    #    plt.hist(np.ravel(test_data),50,normed=True)
    #    plt.plot(x,Rayleigh(x,popt[0]),'g-',linewidth=3)

    #    chi_sq=np.sum((Rayleigh(x,popt[0])-density)**2/sigma**2)
    #    print "Chisq:", chi_sq, chi_sq/(len(hist)-1)
    #    plt.show()


if __name__ == "__main__":
    noise_map()
