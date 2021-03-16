"""Plot the PSF of a wsclean image."""

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import pyplot as plt

with fits.open("uvdump_-psf.fits") as hdus:
    hdu = hdus[0]
    hdr = hdu.header
    data = hdu.data
    #  data[0, 0, :, :]


wcs = WCS(hdr)

center_pix = 1024
width = 32

#  ra, dec, fr, s = wcs.all_pix2world(1024, 1024, 0, 0, 0)

plt.style.use("seaborn")

ax = plt.subplot(projection=wcs[0, 0, :, :])
im = ax.imshow(
    data[
        0,
        0,
        int(center_pix - width / 2) : int(center_pix + width / 2),
        int(center_pix - width / 2) : int(center_pix + width / 2),
    ],
    origin="lower",
    cmap="viridis",
)
ax.set_title(f"MWA PSF @ {hdr['CRVAL3']} Hz")
ax.grid(True, color="white", ls="dotted", alpha=0.4)
plt.savefig("psf.png")
plt.close()

print(data.shape)

ra_pix = np.arange(data.shape[3])
dec_pix = np.arange(data.shape[2])

ra, _, _, _ = wcs.all_pix2world(ra_pix, 0, 0, 0, 0)
_, dec, _, _ = wcs.all_pix2world(0, dec_pix, 0, 0, 0)

plt.plot(
    dec[int(center_pix - width / 2) : int(center_pix + width / 2)],
    data[0, 0, int(center_pix - width / 2) : int(center_pix + width / 2), center_pix],
    lw=2,
    color="seagreen",
)
plt.title(f"MWA PSF @ {hdr['CRVAL3']} Hz")
plt.xlabel("Declination [deg]")
plt.ylabel("PSF")
plt.tight_layout()
plt.savefig("psf_dec.png")
