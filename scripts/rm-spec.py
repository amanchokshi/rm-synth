import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from matplotlib import pyplot as plt

# Read the POGs tables in to return an astropy Table object
exgal = Table.read("../slurm/iono/POGS-II_ExGal.fits")
df_eg = exgal.to_pandas()
df_eg["catalog_id"] = df_eg["catalog_id"].str.decode("utf-8")
df_eg.set_index("catalog_id", inplace=True)

pogs_ii_024 = df_eg.loc["POGSII-EG-024"]
rm_peak = SkyCoord(pogs_ii_024.ra, pogs_ii_024.dec, unit="deg")

with fits.open("../data/1086351512/rts_imgr/run_i/cubes/rts_imgr_p.phi.dirty.fits") as hdus:
    hdu = hdus[0]
    hdr = hdu.header
    data_p = hdu.data

with fits.open("../data/1086351512/rts_imgr/run_i/cubes/rts_imgr_q.phi.dirty.fits") as hdus:
    hdu = hdus[0]
    #  hdr = hdu.header
    data_q = hdu.data

with fits.open("../data/1086351512/rts_imgr/run_i/cubes/rts_imgr_u.phi.dirty.fits") as hdus:
    hdu = hdus[0]
    hdr = hdu.header
    data_u = hdu.data

wcs = WCS(hdr)
phi, _, _ = wcs.all_pix2world(np.arange(data_p.shape[2]), 0, 0, 0)
_, ras, _ = wcs.all_pix2world(0, np.arange(data_p.shape[1]), 0, 0)
_, _, decs = wcs.all_pix2world(0, 0, np.arange(data_p.shape[0]), 0)


_, rm_x, rm_y = wcs.all_world2pix(0, rm_peak.ra.deg, rm_peak.dec.deg, 0)


#  plt.style.use("seaborn")

#  ax = plt.subplot(projection=wcs[:, :, 0])
#  im = ax.imshow(data[:, :, 411], origin="lower", cmap="viridis",)
#  ax.plot(rm_x, rm_y, "darkorange", marker="o", mfc="none", ms=36, mew=2.0)


#  ax.coords.grid(True, color="white", ls="dotted")
#  ax.coords[0].set_format_unit(u.deg)
#  ax.coords[0].set_auto_axislabel(False)
#  ax.set_xlabel("Right Ascension [deg]")

#  ax.coords[1].set_format_unit(u.deg)
#  ax.coords[1].set_auto_axislabel(False)
#  ax.set_ylabel("Declination [deg]")
#  plt.show()

plt.style.use("seaborn")

plt.plot(phi, data_p[128, 128, :], label="p")
#  plt.plot(phi, data_q[128, 128, :], label="q")
#  plt.plot(phi, data_u[128, 128, :], label="u")
plt.xlim([-50, 50])
plt.xlabel("Faraday Depth")
plt.ylabel("Polarized Flux Density")
plt.legend()
plt.tight_layout()
plt.savefig("rm_spec.png")
