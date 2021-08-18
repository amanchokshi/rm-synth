import sys

import numpy as np
from matplotlib import pyplot as plt

sys.path.append("..")
import rm_theory as rmt


def set_size(width, fraction=1):
    """ Set aesthetic figure dimensions to avoid scaling in latex.

    Parameters
    ----------
    width: float
            Width in pts
    fraction: float
            Fraction of the width which you wish the figure to occupy

    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    # Width of figure
    fig_width_pt = width * fraction

    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    golden_ratio = (5 ** 0.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio

    return fig_width_in, fig_height_in


if __name__ == "__main__":

    # Epic spectral colours
    colours = ["#DA3752",
               "#FCAD61",
               "#66C1A4",
               "#3287BC",
               "#5E4FA1"]

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
    I_0, Q_0, U_0, V_0 = rmt.get_IQUV_complex(
        freqs, rm, ref_I_Jy, ref_V_Jy, SI, frac_pol, ref_chi=0.0, ref_freq=200e6
    )

    fdfs = []
    Is = []
    Qs = []
    Us = []
    Vs = []

    G_Ax = 1 + 0j
    G_Ay = 1 + 0j
    G_Bx = 1 + 0j
    G_By = 1 + 0j

    l_Ax = 0 + 0j
    l_Ay = 0 + 0j
    l_Bx = 0 + 0j
    l_By = 0 + 0j

    I, Q, U, V = rmt.rm_leakage(
        I_0, Q_0, U_0, V_0, G_Ax, G_Ay, G_Bx, G_By, l_Ax, l_Ay, l_Bx, l_By
    )

    # Convert stokes to instrumental pols
    XX, XY, YX, YY = rmt.stokes_instru(I, Q, U, V)

    # Determine FDF, RMSF
    fdf, rmsf, phi = rmt.rm_synth(freqs, Q, U, phi_lim=200, dphi=0.1)
    Is.append(I)
    Qs.append(Q)
    Us.append(U)
    Vs.append(V)

    fdfs.append(fdf)

    a = 0.9
    #  p = np.deg2rad(40)

    #  for p in np.arange(0, 225, 45):
    for p in [0, 90, 180]:
        #  for a in np.linspace(0.8, 1.0, 11)[::-1]:

        G = a * np.exp(1j * np.deg2rad(p))

        # Apply some Hamaker leakage
        G_Ax = G
        G_Ay = 1 + 0j
        G_Bx = G
        G_By = 1 + 0j

        l_Ax = 0 + 0j
        l_Ay = 0 + 0j
        l_Bx = 0 + 0j
        l_By = 0 + 0j

        I, Q, U, V = rmt.rm_leakage(
            I_0, Q_0, U_0, V_0, G_Ax, G_Ay, G_Bx, G_By, l_Ax, l_Ay, l_Bx, l_By
        )

        # Convert stokes to instrumental pols
        XX, XY, YX, YY = rmt.stokes_instru(I, Q, U, V)

        # Determine FDF, RMSF
        fdf, rmsf, phi = rmt.rm_synth(freqs, Q, U, phi_lim=200, dphi=0.1)
        fdfs.append(fdf)
        Is.append(I)
        Qs.append(Q)
        Us.append(U)
        Vs.append(V)

    #############################################

    ls = ["solid", "dotted", "dashed", "dashdot"]

    plt.style.use("seaborn")
    plt.rcParams.update(
        {
            #  "font.size": 15,
            "text.usetex": True,
            "font.family": "serif",
            #  "font.serif": "Times New Roman",
            # Use 10pt font in plots, to match 10pt font in document
            "axes.labelsize": 8,
            "axes.titlesize": 9,
            "font.size": 8,
            # Make the legend/label fonts a little smaller
            "legend.fontsize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
        }
    )

    fig = plt.figure(figsize=(7.6, 9.4))
    #  plt.subplots_adjust(wspace=0.35, hspace=0.25)

    #  ax1 = fig.add_subplot(2, 1, 1)
    ax1 = fig.add_axes([0.08, 0.6, 0.9, 0.38])

    ax1.plot(phi, np.abs(fdfs[0]), color="#222222", lw=1, ls=ls[0], label="Simulation")
    ax1.plot(
        phi,
        np.abs(fdfs[1]),
        color=colours[4],
        lw=1.4,
        ls=ls[0],
        label=r"A:0.9, $\phi:0^\circ$",
    )
    ax1.plot(
        phi,
        np.abs(fdfs[2]),
        color=colours[0],
        lw=1.4,
        ls=ls[0],
        label=r"A:0.9, $\phi:90^\circ$",
    )
    ax1.plot(
        phi,
        np.abs(fdfs[3]),
        color=colours[1],
        lw=1.4,
        ls=ls[0],
        label=r"A:0.9, $\phi:180^\circ$",
        zorder=1,
    )

    ax1.set_xlabel(r"Faraday Depth [rad/m$^2$]")
    ax1.set_ylabel("Polarised Flux Density [Jy/PSF/RMSF]")
    ax1.set_xlim([-50, 50])
    leg = ax1.legend(loc="upper right", frameon=True)
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    #  ax2 = fig.add_subplot(4, 2, 5)
    ax2 = fig.add_axes([0.08, 0.31, 0.41, 0.23])
    ax2.plot(
        freqs / 1e6, Is[0].real, color="#222222", lw=1, ls=ls[0], label="Simulation"
    )
    ax2.plot(
        freqs / 1e6,
        Is[1].real,
        color=colours[4],
        lw=1.4,
        ls=ls[0],
        label=r"A:0.9, $\phi:0^\circ$",
    )
    ax2.plot(
        freqs / 1e6,
        Is[2].real,
        color=colours[0],
        lw=1.4,
        ls=ls[0],
        label=r"A:0.9, $\phi:90^\circ$",
    )
    ax2.plot(
        freqs / 1e6,
        Is[3].real,
        color=colours[1],
        lw=1.4,
        ls=ls[0],
        label=r"A:0.9, $\phi:180^\circ$",
        zorder=1,
    )

    ax2.set_ylabel("Stokes I Flux [Jy]")
    ax2.set_xticklabels([])

    #  ax3 = fig.add_subplot(4, 2, 6)
    ax3 = fig.add_axes([0.57, 0.31, 0.41, 0.23])
    ax3.plot(
        freqs / 1e6, Qs[0].real, color="#222222", lw=1, ls=ls[0], label="Simulation"
    )
    ax3.plot(
        freqs / 1e6,
        Qs[1].real,
        color=colours[4],
        lw=1.4,
        ls=ls[0],
        label=r"A:0.9, $\phi:0^\circ$",
    )
    ax3.plot(
        freqs / 1e6,
        Qs[2].real,
        color=colours[0],
        lw=1.4,
        ls=ls[0],
        label=r"A:0.9, $\phi:90^\circ$",
    )
    ax3.plot(
        freqs / 1e6,
        Qs[3].real,
        color=colours[1],
        lw=1.4,
        ls=ls[0],
        label=r"A:0.9, $\phi:180^\circ$",
        zorder=1,
    )

    ax3.set_ylabel("Stokes Q Flux [Jy]")
    ax3.set_xticklabels([])

    #  ax4 = fig.add_subplot(4, 2, 7)
    ax4 = fig.add_axes([0.08, 0.06, 0.41, 0.23])
    ax4.plot(
        freqs / 1e6, Us[0].real, color="#222222", lw=1, ls=ls[0], label="Simulation"
    )
    ax4.plot(
        freqs / 1e6,
        Us[1].real,
        color=colours[4],
        lw=1.4,
        ls=ls[0],
        label=r"A:0.9, $\phi:0^\circ$",
    )
    ax4.plot(
        freqs / 1e6,
        Us[2].real,
        color=colours[0],
        lw=1.4,
        ls=ls[0],
        label=r"A:0.9, $\phi:90^\circ$",
    )
    ax4.plot(
        freqs / 1e6,
        Us[3].real,
        color=colours[1],
        lw=1.4,
        ls=ls[0],
        label=r"A:0.9, $\phi:180^\circ$",
        zorder=1,
    )

    ax4.set_xlabel(r"Frequency [MHz]")
    ax4.set_ylabel("Stokes U Flux [Jy]")

    #  ax5 = fig.add_subplot(4, 2, 8)
    ax5 = fig.add_axes([0.57, 0.06, 0.41, 0.23])
    ax5.plot(
        freqs / 1e6, Vs[0].real, color="#222222", lw=1, ls=ls[0], label="Simulation"
    )
    ax5.plot(
        freqs / 1e6,
        Vs[1].real,
        color=colours[4],
        lw=1.4,
        ls=ls[0],
        label=r"A:0.9, $\phi:0^\circ$",
    )
    ax5.plot(
        freqs / 1e6,
        Vs[2].real,
        color=colours[0],
        lw=1.4,
        ls=ls[0],
        label=r"A:0.9, $\phi:90^\circ$",
    )
    ax5.plot(
        freqs / 1e6,
        Vs[3].real,
        color=colours[1],
        lw=1.4,
        ls=ls[0],
        label=r"A:0.9, $\phi:180^\circ$",
        zorder=1,
    )

    ax5.set_xlabel(r"Frequency [MHz]")
    ax5.set_ylabel("Stokes V Flux [Jy]")

    #  fig.tight_layout()
    #  plt.show()
    plt.savefig("./plots/rm_theory_leakage.pdf", bbox_inches="tight", dpi=300)
    plt.close()
