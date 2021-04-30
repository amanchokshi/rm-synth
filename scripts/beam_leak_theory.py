"""Plot unpolarised all-sky MWA FEE beam response using hyperbeam."""

import concurrent.futures
from itertools import repeat

import healpy as hp
import mwa_hyperbeam
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from tqdm import tqdm

import healpix as hpx
import rm_theory as rmt


def makeUnpolInstrumentalResponse(j1, j2):
    """Convert Jones matracies to unpolarised beam responses.

       Form the visibility matrix in instrumental response from two Jones
       matrices assuming unpolarised sources (hence the brightness matrix is
       the identity matrix)

       Input: j1,j2: Jones matrices of dimension[za][az][2][2]
       Returns: [za][az][[xx,xy],[yx,yy]] where "X" and "Y" are defined by the receptors
       of the Dipole object used in the ApertureArray. Hence to get "XX", you want
       result[za][az][0][0] and for "YY" you want result[za][az][1][1]

       Modified from:
       https://github.com/MWATelescope/mwa_pb/blob/696c2835f44de510da2d5d5bcd3c15223bbe6d7b/mwa_pb/beam_tools.py

       output: [XX, XY, YX, YY]
    """
    result = np.empty_like(j1)

    result[:, 0] = j1[:, 0] * j2[:, 0].conjugate() + j1[:, 1] * j2[:, 1].conjugate()
    result[:, 3] = j1[:, 2] * j2[:, 2].conjugate() + j1[:, 3] * j2[:, 3].conjugate()
    result[:, 1] = j1[:, 0] * j2[:, 2].conjugate() + j1[:, 1] * j2[:, 3].conjugate()
    result[:, 2] = j1[:, 2] * j2[:, 0].conjugate() + j1[:, 3] * j2[:, 1].conjugate()
    return result


if __name__ == "__main__":

    # We can make a new beam object with a path to the HDF5 file specified by MWA_BEAM_FILE.
    beam = mwa_hyperbeam.FEEBeam()

    # Create a grid across the sky at which to evaluate the beam
    #  nside = 501
    nside = 32

    #  x_range = np.linspace(-np.pi / 2, np.pi / 2, nside)
    #  y_range = np.linspace(-np.pi / 2, np.pi / 2, nside)

    #  x_mesh, y_mesh = np.meshgrid(x_range, y_range)

    #  az_grid = np.pi - np.arctan2(x_mesh, y_mesh)
    #  za_grid = np.sqrt(x_mesh ** 2 + y_mesh ** 2)

    #  below_horizon = np.where(np.sqrt(x_mesh ** 2 + y_mesh ** 2) > np.pi / 2)

    #  az_grid[below_horizon] = np.nan
    #  za_grid[below_horizon] = np.nan

    #  az_f = az_grid.flatten()
    #  za_f = za_grid.flatten()

    za, az = hpx.healpix_za_az(nside=nside)

    freq = 200e6
    delays = [0] * 16
    amps_16 = [1.0] * 16
    norm_to_zenith = True

    # beam.calc_jones is also available, but that makes a single Jones matrix at a
    # time, so one would need to iterate over az and za. calc_jones_array is done in
    # parallel with Rust (so it's fast).
    jones_16 = beam.calc_jones_array(az, za, freq, delays, amps_16, norm_to_zenith)

    # Loop over all dead dipoles
    amps_f = np.array([1] + [0] * 15)
    for j in range(16):

        print(f"Dipole Flagged : [{j}]")
        amps_flag = np.roll(amps_f, j)

        jones_15 = beam.calc_jones_array(
            az, za, freq, delays, amps_flag, norm_to_zenith
        )

        # Difference b/w jones matracies added to a unit matrix which
        # represents the perfect response of the instrument
        jones_err = jones_16 - jones_15 + np.array([1.0, 0.0, 0.0, 1.0])

        # MWA constants
        low_freq = 160e6
        high_freq = 230e6
        fine_channel = 40e3

        # RM Source Constants
        SI = -0.7
        rm = 20
        frac_pol = 0.20
        ref_I_Jy = 7
        ref_V_Jy = 1

        # Arrays of freqencies in Hz
        freqs = np.arange(low_freq, high_freq + fine_channel, fine_channel)

        # Get stokes parameters as a function of frequency
        I, Q, U, V = rmt.get_IQUV_complex(
            freqs, rm, ref_I_Jy, ref_V_Jy, SI, frac_pol, ref_chi=0.0, ref_freq=200e6
        )

        beam_loop = True
        if beam_loop:

            fractional_leakage = []

            #  for i, j_err in enumerate(jones_err):
            for i in tqdm(range(len(jones_err))):

                j_err = jones_err[i]

                #  if nan in j_err:
                if np.isnan(np.min(j_err)):
                    fractional_leakage.append(np.nan)
                else:
                    # Apply Hamaker leakage
                    G_Ax = j_err[0]
                    G_Ay = 1 + 0j
                    G_Bx = j_err[3]
                    G_By = 1 + 0j

                    l_Ax = j_err[1]
                    l_Ay = 0 + 0j
                    l_Bx = j_err[2]
                    l_By = 0 + 0j

                    I, Q, U, V = rmt.rm_leakage(
                        I, Q, U, V, G_Ax, G_Ay, G_Bx, G_By, l_Ax, l_Ay, l_Bx, l_By
                    )

                    # Determine FDF, RMSF
                    fdf, rmsf, phi = rmt.rm_synth(freqs, Q, U, phi_lim=200, dphi=0.1)

                    # Find peak of leakage and FDF
                    leak_flux = np.amax(np.abs(fdf)[np.where((phi > -1) & (phi < 1))])
                    pos_flux = np.amax(
                        np.abs(fdf)[np.where((phi > rm - 1) & (phi < rm + 1))]
                    )
                    neg_flux = np.amax(
                        np.abs(fdf)[np.where((phi > -1 * rm - 1) & (phi < -1 * rm + 1))]
                    )

                    if neg_flux > pos_flux:
                        p = -20
                        rm_flux = -1 * neg_flux
                    else:
                        p = 20
                        rm_flux = pos_flux

                    frac_leak = leak_flux / rm_flux

                    fractional_leakage.append(frac_leak)

                    #  plt.style.use("seaborn")
                    #  plt.plot(phi, np.abs(fdf))
                    #  plt.xlim([-50, 50])
                    #  plt.title(
                    #  rf"FDF : $\phi$=20 rad m$^{{-2}}$, Frac Leakage : {frac_leak:.2f}"
                    #  )
                    #  plt.xlabel("Faraday Depth [rad/m$^2$]")
                    #  plt.ylabel("Polarized Flux Density [Jy/PSF/RMSF]")
                    #  plt.scatter(
                    #  [0, p],
                    #  [leak_flux, rm_flux],
                    #  marker="o",
                    #  color="midnightblue",
                    #  facecolor="none",
                    #  linewidth=2.1,
                    #  zorder=7,
                    #  )
                    #  plt.tight_layout()
                    #  plt.savefig(f"./rm_leak/{i}.png")
                    #  plt.close()

            fractional_leakage = np.asarray(fractional_leakage)
            np.save(f"./rm_leak/beam_leak_{j}.npy", fractional_leakage)

            npix = hp.nside2npix(nside=nside)
            beam_rm_leak = np.zeros(npix)
            above_horizon = range(int(npix / 2))
            beam_rm_leak[above_horizon] = fractional_leakage

            plt.style.use("seaborn")
            fig = plt.figure(figsize=(9, 9))
            fig.suptitle(
                f"Fractional Polarization Leakage - Flagged Dipole: [{j}]",
                fontsize=16,
                y=0.92,
            )
            hpx.plot_healpix(
                data_map=beam_rm_leak, sub=(1, 1, 1), cmap="Spectral",
            )
            plt.savefig(f"./rm_leak/beam_leak_{j}.png", bbox_inches="tight")

            #  plt.style.use("seaborn")
            #  fractional_leakage = fractional_leakage.reshape(az_grid.shape)
            #  fig, ax = plt.subplots()
            #  ax.grid(True, color="white", alpha=0.7, ls="dotted")
            #  im = ax.imshow(fractional_leakage, cmap="viridis")
            #  cbar = ax.figure.colorbar(im, label="Fractional Leakage")
            #  plt.tight_layout()
            #  #  plt.savefig("beam_leak.png")
            #  plt.show()

    #############################################################################
    #                                                                           #
    #                       Healpix XX, YY Beam Maps                            #
    #                                                                           #
    #############################################################################

    #  za, az = hpx.healpix_za_az(nside=32)
    #  jn = beam.calc_jones_array(az, za, freq, delays, amps_16, norm_to_zenith)
    #  unpol_beam = makeUnpolInstrumentalResponse(jn, jn)

    #  XX = unpol_beam[:, 0]
    #  YY = unpol_beam[:, 3]

    #  npix = hp.nside2npix(nside=32)
    #  beam_response_XX = np.zeros(npix)
    #  beam_response_YY = np.zeros(npix)
    #  above_horizon = range(int(npix / 2))
    #  beam_response_XX[above_horizon] = XX
    #  beam_response_YY[above_horizon] = YY

    #  plt.style.use("seaborn")
    #  fig = plt.figure(figsize=(6, 6))
    #  fig.suptitle(f"MWA HYPERBEAM FEE MAP XX", fontsize=16, y=0.92)
    #  hpx.plot_healpix(
    #  data_map=10 * np.log10(beam_response_XX),
    #  sub=(1, 1, 1),
    #  cmap="viridis",
    #  vmin=-50,
    #  vmax=0,
    #  )
    #  plt.savefig("XX.png", bbox_inches="tight")

    #  plt.style.use("seaborn")
    #  fig = plt.figure(figsize=(6, 6))
    #  fig.suptitle(f"MWA HYPERBEAM FEE MAP YY", fontsize=16, y=0.92)
    #  hpx.plot_healpix(
    #  data_map=10 * np.log10(beam_response_YY),
    #  sub=(1, 1, 1),
    #  cmap="viridis",
    #  vmin=-50,
    #  vmax=0,
    #  )
    #  plt.savefig("YY.png", bbox_inches="tight")

    #############################################################################

    #  j_err = jones_err[272]

    #  # Apply Hamaker leakage
    #  G_Ax = j_err[0]
    #  G_Ay = 1 + 0j
    #  G_Bx = j_err[3]
    #  G_By = 1 + 0j

    #  l_Ax = j_err[1]
    #  l_Ay = 0 + 0j
    #  l_Bx = j_err[2]
    #  l_By = 0 + 0j

    #  I, Q, U, V = rmt.rm_leakage(
    #  I, Q, U, V, G_Ax, G_Ay, G_Bx, G_By, l_Ax, l_Ay, l_Bx, l_By
    #  )

    #  # Determine FDF, RMSF
    #  fdf, rmsf, phi = rmt.rm_synth(freqs, Q, U, phi_lim=200, dphi=0.1)

    #  # Find peak of leakage and FDF
    #  leak_flux = np.amax(np.abs(fdf)[np.where((phi > -1) & (phi < 1))])
    #  pos_flux = np.amax(np.abs(fdf)[np.where((phi > rm - 1) & (phi < rm + 1))])
    #  neg_flux = np.amax(
    #  np.abs(fdf)[np.where((phi > -1 * rm - 1) & (phi < -1 * rm + 1))]
    #  )

    #  if neg_flux > pos_flux:
    #  p = -20
    #  rm_flux = -1 * neg_flux
    #  else:
    #  p = 20
    #  rm_flux = pos_flux

    #  frac_leak = leak_flux / rm_flux

    #  plt.style.use("seaborn")
    #  plt.plot(phi, np.abs(fdf))
    #  plt.xlim([-50, 50])
    #  plt.title(
    #  rf"FDF : $\phi$=20 rad m$^{{-2}}$, Frac Leakage : {frac_leak:.2f}"
    #  )
    #  plt.xlabel("Faraday Depth [rad/m$^2$]")
    #  plt.ylabel("Polarized Flux Density [Jy/PSF/RMSF]")
    #  plt.scatter(
    #  [0, p],
    #  [leak_flux, rm_flux],
    #  marker="o",
    #  color="midnightblue",
    #  facecolor="none",
    #  linewidth=2.1,
    #  zorder=7,
    #  )

    #  plt.tight_layout()
    #  plt.show()

    ####################################################################

    #  unpol_beam = makeUnpolInstrumentalResponse(jones_err, jones_err)

    #  XX = unpol_beam[:, 0]
    #  YY = unpol_beam[:, 3]

    #  beam_weights = (XX + YY) / 2.0

    #  beam_weights = beam_weights.reshape(az_grid.shape)

    #  #  plt.style.use("seaborn")
    #  #  plt.imshow(beam_weights.real, norm=LogNorm(vmin=0.0001, vmax=1))
    #  #  plt.tight_layout()
    #  #  plt.show()

    #  levels = [0.001, 0.01, 0.1, 0.3, 0.5, 0.9]

    #  fig, ax = plt.subplots()
    #  im = ax.imshow(beam_weights.real, norm=LogNorm(vmin=0.0001, vmax=1))
    #  CS = ax.contour(beam_weights.real, levels, colors="white", linewidths=0.7)
    #  ax.clabel(CS, inline=1, fontsize=5)
    #  plt.tight_layout()
    #  plt.show()

    #####################################################################
