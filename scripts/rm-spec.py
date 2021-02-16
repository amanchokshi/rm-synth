"""
Plot the RM Specta of a POGS source.
"""

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from matplotlib import pyplot as plt


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
    Astropy.coordinates.SkyCoord object
        Ra, Dec of POGS object in radians
    """

    exgal = Table.read(pogs_path)
    df_eg = exgal.to_pandas()
    df_eg["catalog_id"] = df_eg["catalog_id"].str.decode("utf-8")
    df_eg.set_index("catalog_id", inplace=True)

    pogs_source = df_eg.loc[f"{pogs_obj}"]
    pogs_pos = SkyCoord(pogs_source.ra, pogs_source.dec, unit="deg")

    return pogs_pos


#  # Read the POGs tables in to return an astropy Table object
#  exgal = Table.read("../slurm/iono/POGS-II_ExGal.fits")
#  df_eg = exgal.to_pandas()
#  df_eg["catalog_id"] = df_eg["catalog_id"].str.decode("utf-8")
#  df_eg.set_index("catalog_id", inplace=True)

#  pogs = df_eg.loc["POGSII-EG-321"]
#  rm_peak = SkyCoord(pogs.ra, pogs.dec, unit="deg")

rm_peak = pogs_obj_loc("POGSII-EG-321", "../slurm/iono/POGS-II_ExGal.fits")


##########################
# data == [dec, ra, phi] #
##########################

with fits.open("../data/rts_imgr_p.phi.dirty.fits") as hdus:
    hdu = hdus[0]
    hdr = hdu.header
    data_p = hdu.data

print(repr(hdr))

#  with fits.open("../data/1086351512/rts_imgr/run_i/cubes/rts_imgr_q.phi.dirty.fits") as hdus:
#  hdu = hdus[0]
#  #  hdr = hdu.header
#  data_q = hdu.data

#  with fits.open("../data/1086351512/rts_imgr/run_i/cubes/rts_imgr_u.phi.dirty.fits") as hdus:
#  hdu = hdus[0]
#  hdr = hdu.header
#  data_u = hdu.data

wcs = WCS(hdr)
phi, _, _ = wcs.all_pix2world(np.arange(data_p.shape[2]), 0, 0, 0)
_, ras, _ = wcs.all_pix2world(0, np.arange(data_p.shape[1]), 0, 0)
_, _, decs = wcs.all_pix2world(0, 0, np.arange(data_p.shape[0]), 0)


_, ra_x, dec_y = wcs.all_world2pix(0, rm_peak.ra.deg, rm_peak.dec.deg, 0)

ra_x = np.round(ra_x).astype(int)
dec_y = np.round(dec_y).astype(int)

spec = data_p[np.round(dec_y).astype(int), np.round(ra_x).astype(int), :]
spec_peak = np.amax(spec)
peak_idx = np.where(spec == spec_peak)[0][0]

print(f"Phi Max: {phi[peak_idx]}, RM: {spec_peak:.3f}")

plt.style.use("seaborn")

ax = plt.subplot(projection=wcs[:, :, 355])
im = ax.imshow(data_p[:, :, 355], origin="lower", cmap="viridis",)
ax.plot(
    np.round(ra_x).astype(int),
    np.round(dec_y).astype(int),
    "darkorange",
    marker="o",
    mfc="none",
    ms=28,
    mew=2.0,
)

ax.set_xlim(ra_x - 50, ra_x + 50)
ax.set_ylim(dec_y - 50, dec_y + 50)


ax.coords.grid(True, color="white", ls="dotted")
ax.coords[0].set_format_unit(u.deg)
ax.coords[0].set_auto_axislabel(False)
ax.set_xlabel("Right Ascension [deg]")

ax.coords[1].set_format_unit(u.deg)
ax.coords[1].set_auto_axislabel(False)
ax.set_ylabel("Declination [deg]")

plt.show()

plt.style.use("seaborn")

plt.plot(
    phi, data_p[np.round(dec_y).astype(int), np.round(ra_x).astype(int), :], label="p"
)
#  plt.plot(phi, data_q[128, 128, :], label="q")
#  plt.plot(phi, data_u[128, 128, :], label="u")
plt.xlabel("Faraday Depth")
plt.ylabel("Polarized Flux Density")
plt.legend()
plt.tight_layout()
#  plt.savefig("rm_spec.png")
plt.show()
