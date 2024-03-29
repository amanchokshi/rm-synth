import numpy as np
from matplotlib import pyplot as plt
from scipy import constants as const


def spectal_index(freqs, SI, ref_flux, ref_freq):
    """Determine flux at frequencies based on SI and ref_freq.

    Parameters
    ----------
    freqs : numpy.array / float
        An array of freqencies in Hz at which flux is to be determined
    SI : float
        Spectral index of source
    ref_flux : float
        Reference flux in Jy
    ref_freq : float
        Reference frequency in HZ

    Returns
    -------
    numpy.array / float
        flux_freqs : An array of fluxes in Jy, for every frequency
    """

    flux_freqs = ref_flux * ((freqs / ref_freq) ** SI)

    return flux_freqs


def pol_angle_lamba(rm, lambdas, ref_chi):
    """Polarization angle as a function of frequency.

    Parameters
    ----------
    rm : float
        Rotation measure of source in [rad m^-2]
    lambdas : numpy.array / float
        Wavelengths in metres
    ref_chi : float
        Reference polarization angle

    Returns
    -------
    numpy.array / float
        Polarization angle at wavelengths
    """

    return ref_chi + rm * lambdas ** 2


def get_IQUV_complex(
    freqs, rm, ref_I_Jy, ref_V_Jy, SI, frac_pol, ref_chi=0.0, ref_freq=200e6
):
    """
    Get I, Q, U stokes parameters as a function of freqency.

    Parameters
    ----------
    freqs : numpy.array / float
        An array of freqencies in Hz at which flux is to be determined
    rm : float
        Rotation measure of source in [rad m^-2]
    ref_I_Jy : float
        Reference flux in Jy of stokes I
    SI : float
        Spectral index of source
    frac_pol : float
        Polarization fraction
    ref_chi : float
        Reference polarization angle, default: 0.0
    ref_freq : float
        Reference frequency in HZ, default : 200e6

    Returns
    -------
    I, Q, U complex stokes parameters as a function of frequency

    Note
    ----
    Modified from Jack Line's beam test script
    """
    lambdas = const.c / freqs

    pol_ang = pol_angle_lamba(rm, lambdas, ref_chi)

    I = spectal_index(freqs, SI, ref_I_Jy, ref_freq)

    numer = frac_pol * I * np.exp(2j * pol_ang)
    denom = 1 + 1j * np.tan(2 * pol_ang)

    Q = numer / denom

    U = Q * np.tan(2 * pol_ang)

    V = spectal_index(freqs, SI, ref_V_Jy, ref_freq=200e6)

    return I, Q, U, V


def stokes_instru(I, Q, U, V):
    """convert stokes parameters to instrumental frame.

    parameters
    ----------
    i, q, u, v : numpy.array
        complex stokes parameters

    returns
    -------
    xx, yy, xy, yx : numpy.array
        instrumental polarized measurement
    """

    XX = I + Q
    YY = I - Q
    XY = U + 1j * V
    YX = U - 1j * V

    return XX, XY, YX, YY


def instru_stokes(XX, XY, YX, YY):
    """convert instrumental measurement to stokes parameters.

    parameters
    ----------
    XX, YY, XY, YX : numpy.array
        instrumental polarized measurement

    returns
    -------
    I, Q, U, V : numpy.array
        complex stokes parameters
    """

    I = (XX + YY) / 2
    Q = (XX - YY) / 2
    U = (XY + YX) / 2
    V = (1j * (XY + YX)) / 2

    return I, Q, U, V


def rm_leakage(I, Q, U, V, G_Ax, G_Ay, G_Bx, G_By, l_Ax, l_Ay, l_Bx, l_By):
    """Hamaker based leakage."""

    # Complex conjugate of Jones matrix B
    G_Bx_c = np.conj(G_Bx)
    G_By_c = np.conj(G_By)
    l_Bx_c = np.conj(l_Bx)
    l_By_c = np.conj(l_By)

    # Stokes parameters with leakage
    # Eqns 4, 5, 6, 7 in the rm_theory.ipynb
    I_m = 0.5 * (
        (G_Ax * G_Bx_c + l_Ax * l_Bx_c + l_Ay * l_By_c + G_Ay * G_By_c) * I
        + (G_Ax * G_Bx_c - l_Ax * l_Bx_c + l_Ay * l_By_c - G_Ay * G_By_c) * Q
        + (G_Ax * l_Bx_c + l_Ax * G_Bx_c - l_Ay * G_By_c - G_Ay * l_By_c) * U
        + (
            1j * G_Ax * l_Bx_c
            - 1j * l_Ax * G_Bx_c
            - 1j * l_Ay * G_By_c
            + 1j * G_Ay * l_By_c
        )
        * V
    )

    Q_m = 0.5 * (
        (G_Ax * G_Bx_c + l_Ax * l_Bx_c - l_Ay * l_By_c - G_Ay * G_By_c) * I
        + (G_Ax * G_Bx_c - l_Ax * l_Bx_c - l_Ay * l_By_c + G_Ay * G_By_c) * Q
        + (G_Ax * l_Bx_c + l_Ax * G_Bx_c + l_Ay * G_By_c + G_Ay * l_By_c) * U
        + (
            1j * G_Ax * l_Bx_c
            - 1j * l_Ax * G_Bx_c
            + 1j * l_Ay * G_By_c
            - 1j * G_Ay * l_By_c
        )
        * V
    )

    U_m = 0.5 * (
        (-G_Ax * l_By_c + l_Ax * G_By_c - l_Ay * G_Bx_c + G_Ay * l_Bx_c) * I
        + (-G_Ax * l_By_c - l_Ax * G_By_c - l_Ay * G_Bx_c - G_Ay * l_Bx_c) * Q
        + (G_Ax * G_By_c - l_Ax * l_By_c - l_Ay * l_Bx_c + G_Ay * G_Bx_c) * U
        + (
            1j * G_Ax * G_By_c
            + 1j * l_Ax * l_By_c
            - 1j * l_Ay * l_Bx_c
            - 1j * G_Ay * G_Bx_c
        )
        * V
    )
    V_m = -0.5j * (
        (-G_Ax * l_By_c + l_Ax * G_By_c + l_Ay * G_Bx_c - G_Ay * l_Bx_c) * I
        + (-G_Ax * l_By_c - l_Ax * G_By_c + l_Ay * G_Bx_c + G_Ay * l_Bx_c) * Q
        + (G_Ax * G_By_c - l_Ax * l_By_c + l_Ay * l_Bx_c - G_Ay * G_Bx_c) * U
        + (
            1j * G_Ax * G_By_c
            + 1j * l_Ax * l_By_c
            + 1j * l_Ay * l_Bx_c
            + 1j * G_Ay * G_Bx_c
        )
        * V
    )

    return I_m, Q_m, U_m, V_m


def rm_synth(freqs, Q, U, phi_lim=200, dphi=0.5):
    """Do RM Synthesis on stokes Q & U vectors.

    Parameters
    ----------
    freqs : numpy.array / float
        An array of freqencies in Hz at which flux is to be determined
    Q, U : numpy.array
        Linear stokes vectors
    phi_lim : float
        Faraday depth limit, default : +-200
    dphi : float
        Faraday depth resolution, default : 0.5

    Returns
    -------
    numpy.array
        FDF (Faraday Dispersion Function)
        RMSF (Rotation Measure Spread Function)
        Phi (Array of faraday depths)
    """

    # Wavelengths
    lambdas = const.c / freqs

    # Uniform weights
    weights = np.ones(Q.shape)

    # Eqn 38 (B&dB 2005)
    K = 1 / np.sum(weights)

    # Eqn 32 (B&dB 2005) - weighted mean on lambda^2
    lam0sq = K * np.sum(weights * lambdas ** 2)

    # Phi array - range of faraday depths
    phi_arr = np.arange(-phi_lim, phi_lim + dphi, dphi)

    # Complex linear polarization
    P = Q + 1j * U

    FDF = K * np.sum(
        P * np.exp(np.outer(-2.0j * phi_arr, (lambdas ** 2 - lam0sq))), axis=1
    )

    RMSF = K * np.sum(
        weights * np.exp(np.outer(-2.0j * phi_arr, (lambdas ** 2 - lam0sq))), axis=1
    )

    return FDF, RMSF, phi_arr


def rm_clean(FDF, RMSF, phi_arr, clean_cut=2):
    """RMCLEAN - Heald et al. 2009"""

    # STD and Median of ||FDF||
    noise = np.std(np.abs(FDF))
    noise_floor = np.median(np.abs(FDF))

    # Will perform rm clean upto clean_cut*noise above the noise_floor
    clean_floor = clean_cut * noise + noise_floor

    print(noise, noise_floor, clean_floor)


def plot_rm_grid(freqs, I, Q, U, V, XX, XY, YX, YY, rmsf, fdf, phi):

    # Plot stokes vectors, FDF, RMSF
    plt.style.use("seaborn")
    plt.figure(figsize=(11, 7))

    # Plot stokes vectors
    ax1 = plt.subplot(221)
    colors = plt.cm.Spectral([0.01, 0.14, 0.86, 0.99])
    ax1.set_ylabel("Flux [Jy]")
    ax1.set_title("Stoke Fluxes vs Frequency")
    stokes = [I, Q, U, V]
    for i, st in enumerate(["I", "Q", "U", "V"]):
        ax1.plot(freqs / 1e6, np.real(stokes[i]), color=colors[i], label=st)

    leg = ax1.legend(frameon=True, markerscale=1, handlelength=1, loc="upper right")
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    # Plot RMSF
    ax2 = plt.subplot(222)
    ax2.set_xlim([-20, 20])
    ax2.plot(phi, np.abs(rmsf), label=r"$ \vert R \vert $", zorder=3)
    ax2.plot(phi, np.real(rmsf), label=r"$ real(R) $")
    ax2.plot(phi, np.imag(rmsf), label=r"$ imag(R) $")
    ax2.set_title("RMSF [-20, 20]")
    ax2.set_ylabel("RMSF")

    leg = ax2.legend(frameon=True, markerscale=1, handlelength=1, loc="upper right")
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    # Plot Instrumental Polarizations
    ax3 = plt.subplot(223)
    colors = plt.cm.Spectral([0.01, 0.14, 0.86, 0.99])
    instru = [XX, YY, XY, YX]
    for i, inst in enumerate(["XX", "YY", "XY", "YX"]):
        ax3.plot(freqs / 1e6, np.real(instru[i]), color=colors[i], label=inst)

    ax3.set_xlabel("Frequency [MHz]")
    ax3.set_ylabel("Flux [Jy]")
    ax3.set_title("Instrumental Pol Fluxes vs Frequency")
    leg = ax3.legend(frameon=True, markerscale=1, handlelength=1, loc="upper right")
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    # Plot FDF
    ax4 = plt.subplot(224)
    ax4.set_xlim([-50, 50])
    ax4.plot(phi, np.abs(fdf), label=r"FDF", zorder=3)
    ax4.set_title(r"FDF : $\phi$=20 rad m$^{-2}$")
    ax4.set_xlabel("Faraday Depth [rad/m$^2$]")
    ax4.set_ylabel("Polarized Flux Density [Jy/PSF/RMSF]")
    leg = ax4.legend(frameon=True, markerscale=1, handlelength=1, loc="upper right")
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":

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
    #  I, Q, U, V = get_IQUV_complex(
    #  freqs, rm, ref_I_Jy, ref_V_Jy, SI, frac_pol, ref_chi=0.0, ref_freq=200e6
    #  )

    #  Q = 0.1 * I + Q
    #  I = 0.9 * I

    # Convert stokes to instrumental pols
    #  XX, XY, YX, YY = stokes_instru(I, Q, U, V)

    # Apply Hamaker leakage

    #  G_Ax = 1 + 0j
    #  G_Ay = 1 + 0j
    #  G_Bx = 1 + 0j
    #  G_By = 1 + 0j

    #  l_Ax = 0 + 0j
    #  l_Ay = 0 + 0j
    #  l_Bx = 0 + 0j
    #  l_By = 0 + 0j

    #  I, Q, U, V = rm_leakage(I, Q, U, V, G_Ax, G_Ay, G_Bx, G_By, l_Ax, l_Ay, l_Bx, l_By)

    # Determine FDF, RMSF
    #  fdf, rmsf, phi = rm_synth(freqs, Q, U, phi_lim=200, dphi=0.1)

    #  rm_clean(fdf, rmsf, phi)

    #  plot_rm_grid(freqs, I, Q, U, V, XX, XY, YX, YY, rmsf, fdf, phi)

    #####################################################################
    #                                                                   #
    #           Save frames of gain amplitude test animation            #
    #                                                                   #
    #####################################################################

    Save = True

    if Save:

        # Get stokes parameters as a function of frequency
        I_0, Q_0, U_0, V_0 = get_IQUV_complex(
            freqs, rm, ref_I_Jy, ref_V_Jy, SI, frac_pol, ref_chi=0.0, ref_freq=200e6
        )

        for fr, a in enumerate(np.linspace(0.8, 1.0, 512)):

            print(f"{fr:04d}: {a:.5f}")

            G = a * np.exp(2j * np.pi)

            # Apply some Hamaker leakage
            G_Ax = G
            G_Ay = 1 + 0j
            G_Bx = G
            G_By = 1 + 0j

            l_Ax = 0 + 0j
            l_Ay = 0 + 0j
            l_Bx = 0 + 0j
            l_By = 0 + 0j

            I, Q, U, V = rm_leakage(
                I_0, Q_0, U_0, V_0, G_Ax, G_Ay, G_Bx, G_By, l_Ax, l_Ay, l_Bx, l_By
            )

            # Convert stokes to instrumental pols
            XX, XY, YX, YY = stokes_instru(I, Q, U, V)

            # Determine FDF, RMSF
            fdf, rmsf, phi = rm_synth(freqs, Q, U, phi_lim=200, dphi=0.1)

            # Plot stokes vectors, FDF, RMSF
            plt.style.use("seaborn")
            fig = plt.figure(figsize=(14, 8))

            # Plot stokes vectors
            ax1 = plt.subplot(221)
            colors = plt.cm.Spectral([0.01, 0.14, 0.86, 0.99])
            ax1.set_ylabel("Flux [Jy]")
            ax1.set_title("Stoke Fluxes vs Frequency")
            stokes = [I, Q, U, V]
            for i, st in enumerate(["I", "Q", "U", "V"]):
                ax1.plot(freqs / 1e6, np.real(stokes[i]), color=colors[i], label=st)

            leg = ax1.legend(
                frameon=True, markerscale=1, handlelength=1, loc="upper right"
            )
            leg.get_frame().set_facecolor("white")
            for le in leg.legendHandles:
                le.set_alpha(1)

            # Plot RMSF
            ax2 = plt.subplot(222)
            ax2.set_xlim([-20, 20])
            ax2.plot(phi, np.abs(rmsf), label=r"$ \vert R \vert $", zorder=3)
            ax2.plot(phi, np.real(rmsf), label=r"$ real(R) $")
            ax2.plot(phi, np.imag(rmsf), label=r"$ imag(R) $")
            ax2.set_title("RMSF [-20, 20]")
            ax2.set_ylabel("RMSF")

            leg = ax2.legend(
                frameon=True, markerscale=1, handlelength=1, loc="upper right"
            )
            leg.get_frame().set_facecolor("white")
            for le in leg.legendHandles:
                le.set_alpha(1)

            # Plot Instrumental Polarizations
            ax3 = plt.subplot(223)
            colors = plt.cm.Spectral([0.01, 0.14, 0.86, 0.99])
            instru = [XX, YY, XY, YX]
            for i, inst in enumerate(["XX", "YY", "XY", "YX"]):
                ax3.plot(freqs / 1e6, np.real(instru[i]), color=colors[i], label=inst)

            ax3.set_xlabel("Frequency [MHz]")
            ax3.set_ylabel("Flux [Jy]")
            ax3.set_title("Instrumental Pol Fluxes vs Frequency")
            leg = ax3.legend(
                frameon=True, markerscale=1, handlelength=1, loc="upper right"
            )
            leg.get_frame().set_facecolor("white")
            for le in leg.legendHandles:
                le.set_alpha(1)

            # Plot FDF
            ax4 = plt.subplot(224)
            ax4.set_xlim([-50, 50])
            ax4.plot(phi, np.abs(fdf), label=r"FDF", zorder=3)
            ax4.set_title(r"FDF : $\phi$=20 rad m$^{-2}$")
            ax4.set_xlabel("Faraday Depth [rad/m$^2$]")
            ax4.set_ylabel("Polarized Flux Density [Jy/PSF/RMSF]")
            leg = ax4.legend(
                frameon=True, markerscale=1, handlelength=1, loc="upper right"
            )
            leg.get_frame().set_facecolor("white")
            for le in leg.legendHandles:
                le.set_alpha(1)

            fig.suptitle(
                f"Amplitude Errors to Gains of X Dipoles : [{a:.5f}]", fontsize=16
            )
            plt.savefig(
                f"../data/rm_theory/gain_amp_2/amp_{fr:04d}_{a:.5f}.png",
                bbox_inches="tight",
                dpi=200,
            )
            plt.close()

    #####################################################################
    #                                                                   #
    #             Save frames of gain phase test animation              #
    #                                                                   #
    #####################################################################

    Save = True

    if Save:

        # Get stokes parameters as a function of frequency
        I_0, Q_0, U_0, V_0 = get_IQUV_complex(
            freqs, rm, ref_I_Jy, ref_V_Jy, SI, frac_pol, ref_chi=0.0, ref_freq=200e6
        )

        for fr, a in enumerate(np.linspace(0, np.pi / 4, 512)):
            #  for fr, a in enumerate(np.linspace(0, 5, 201)):

            deg = np.rad2deg(a)
            print(f"{fr:04d}: {deg:7.3f}")
            #  print(f"{fr:04d}: {a:7.3f}")

            G = 1.0 * np.exp(1j * np.deg2rad(a))

            # Apply some Hamaker leakage
            G_Ax = G
            G_Ay = 1 + 0j
            G_Bx = G
            G_By = 1 + 0j

            l_Ax = 0 + 0j
            l_Ay = 0 + 0j
            l_Bx = 0 + 0j
            l_By = 0 + 0j

            I, Q, U, V = rm_leakage(
                I_0, Q_0, U_0, V_0, G_Ax, G_Ay, G_Bx, G_By, l_Ax, l_Ay, l_Bx, l_By
            )

            # Convert stokes to instrumental pols
            XX, XY, YX, YY = stokes_instru(I, Q, U, V)

            # Determine FDF, RMSF
            fdf, rmsf, phi = rm_synth(freqs, Q, U, phi_lim=200, dphi=0.1)

            # Plot stokes vectors, FDF, RMSF
            plt.style.use("seaborn")

            plt.rcParams.update(
                {
                    #  "font.size": 24,
                    "figure.titlesize": 24,
                    "axes.labelsize": 16,
                    "axes.titlesize": 16,
                    "xtick.labelsize": 14,
                    "ytick.labelsize": 14,
                    "legend.fontsize": 14,
                    "font.family": "serif",
                    "font.serif": "Times New Roman",
                    "text.usetex": True,
                }
            )

            fig = plt.figure(figsize=(14, 8))

            # Plot stokes vectors
            ax1 = plt.subplot(221)
            colors = plt.cm.Spectral([0.01, 0.14, 0.86, 0.99])
            ax1.set_ylabel("Flux [Jy]")
            ax1.set_title("Stoke Fluxes vs Frequency")
            stokes = [I, Q, U, V]
            for i, st in enumerate(["I", "Q", "U", "V"]):
                ax1.plot(freqs / 1e6, np.real(stokes[i]), color=colors[i], label=st)

            leg = ax1.legend(
                frameon=True, markerscale=1, handlelength=1, loc="upper right"
            )
            leg.get_frame().set_facecolor("white")
            for le in leg.legendHandles:
                le.set_alpha(1)

            # Plot RMSF
            ax2 = plt.subplot(222)
            ax2.set_xlim([-20, 20])
            ax2.plot(phi, np.abs(rmsf), label=r"$ \vert R \vert $", zorder=3)
            ax2.plot(phi, np.real(rmsf), label=r"$ real(R) $")
            ax2.plot(phi, np.imag(rmsf), label=r"$ imag(R) $")
            ax2.set_title("RMSF [-20, 20]")
            ax2.set_ylabel("RMSF")

            leg = ax2.legend(
                frameon=True, markerscale=1, handlelength=1, loc="upper right"
            )
            leg.get_frame().set_facecolor("white")
            for le in leg.legendHandles:
                le.set_alpha(1)

            # Plot Instrumental Polarizations
            ax3 = plt.subplot(223)
            colors = plt.cm.Spectral([0.01, 0.14, 0.86, 0.99])
            instru = [XX, YY, XY, YX]
            for i, inst in enumerate(["XX", "YY", "XY", "YX"]):
                ax3.plot(freqs / 1e6, np.real(instru[i]), color=colors[i], label=inst)

            ax3.set_xlabel("Frequency [MHz]")
            ax3.set_ylabel("Flux [Jy]")
            ax3.set_title("Instrumental Pol Fluxes vs Frequency")
            leg = ax3.legend(
                frameon=True, markerscale=1, handlelength=1, loc="upper right"
            )
            leg.get_frame().set_facecolor("white")
            for le in leg.legendHandles:
                le.set_alpha(1)

            # Plot FDF
            ax4 = plt.subplot(224)
            #  ax4 = plt.subplot(111)
            ax4.set_xlim([-50, 50])
            ax4.plot(phi, np.abs(fdf), label=r"FDF", zorder=3)
            # ax4.plot(
            #     phi,
            #     np.abs(fdf),
            #     label=r"FDF",
            #     zorder=3,
            #     linewidth=3.14,
            #     color="seagreen",
            #     alpha=0.9,
            # )
            ax4.set_title(r"FDF : $\phi$=20 rad m$^{-2}$")
            ax4.set_xlabel("Faraday Depth [rad/m$^2$]")
            ax4.set_ylabel("Polarized Flux Density [Jy/PSF/RMSF]")
            leg = ax4.legend(
                frameon=True, markerscale=1, handlelength=1, loc="upper right"
            )
            leg.get_frame().set_facecolor("white")
            for le in leg.legendHandles:
                le.set_alpha(1)

            fig.suptitle(
                f"Phase Errors to Gains of X Dipoles : [{deg:7.3f}$^\circ$]",
                #  f"Phase Errors to Gains of X Dipoles : [{a:7.2f}$^\circ$]",
                fontsize=16,
            )
            plt.savefig(
                f"../data/rm_theory/gain_phase_2/phase_{fr:04d}_{deg:7.3f}.png",
                #  f"../data/phase_sim/tmp2/phase_{fr:04d}.png",
                dpi=200,
                bbox_inches="tight",
                #  transparent=True,
            )
            #  plt.show()
            plt.close()
            #  break

    #####################################################################
    #                                                                   #
    #        Save frames of leakage amplitude test animation            #
    #                                                                   #
    #####################################################################

    Save = True

    if Save:

        # Get stokes parameters as a function of frequency
        I_0, Q_0, U_0, V_0 = get_IQUV_complex(
            freqs, rm, ref_I_Jy, ref_V_Jy, SI, frac_pol, ref_chi=0.0, ref_freq=200e6
        )

        for fr, a in enumerate(np.linspace(0, 0.2, 512)):

            print(f"{fr:04d}: {a:.5f}")

            G = a * np.exp(2j * np.pi)

            # Apply some Hamaker leakage
            G_Ax = 1 + 0j
            G_Ay = 1 + 0j
            G_Bx = 1 + 0j
            G_By = 1 + 0j

            l_Ax = G
            l_Ay = 0 + 0j
            l_Bx = G
            l_By = 0 + 0j

            I, Q, U, V = rm_leakage(
                I_0, Q_0, U_0, V_0, G_Ax, G_Ay, G_Bx, G_By, l_Ax, l_Ay, l_Bx, l_By
            )

            # Convert stokes to instrumental pols
            XX, XY, YX, YY = stokes_instru(I, Q, U, V)

            # Determine FDF, RMSF
            fdf, rmsf, phi = rm_synth(freqs, Q, U, phi_lim=200, dphi=0.1)

            # Plot stokes vectors, FDF, RMSF
            plt.style.use("seaborn")
            fig = plt.figure(figsize=(14, 8))

            # Plot stokes vectors
            ax1 = plt.subplot(221)
            colors = plt.cm.Spectral([0.01, 0.14, 0.86, 0.99])
            ax1.set_ylabel("Flux [Jy]")
            ax1.set_title("Stoke Fluxes vs Frequency")
            stokes = [I, Q, U, V]
            for i, st in enumerate(["I", "Q", "U", "V"]):
                ax1.plot(freqs / 1e6, np.real(stokes[i]), color=colors[i], label=st)

            leg = ax1.legend(
                frameon=True, markerscale=1, handlelength=1, loc="upper right"
            )
            leg.get_frame().set_facecolor("white")
            for le in leg.legendHandles:
                le.set_alpha(1)

            # Plot RMSF
            ax2 = plt.subplot(222)
            ax2.set_xlim([-20, 20])
            ax2.plot(phi, np.abs(rmsf), label=r"$ \vert R \vert $", zorder=3)
            ax2.plot(phi, np.real(rmsf), label=r"$ real(R) $")
            ax2.plot(phi, np.imag(rmsf), label=r"$ imag(R) $")
            ax2.set_title("RMSF [-20, 20]")
            ax2.set_ylabel("RMSF")

            leg = ax2.legend(
                frameon=True, markerscale=1, handlelength=1, loc="upper right"
            )
            leg.get_frame().set_facecolor("white")
            for le in leg.legendHandles:
                le.set_alpha(1)

            # Plot Instrumental Polarizations
            ax3 = plt.subplot(223)
            colors = plt.cm.Spectral([0.01, 0.14, 0.86, 0.99])
            instru = [XX, YY, XY, YX]
            for i, inst in enumerate(["XX", "YY", "XY", "YX"]):
                ax3.plot(freqs / 1e6, np.real(instru[i]), color=colors[i], label=inst)

            ax3.set_xlabel("Frequency [MHz]")
            ax3.set_ylabel("Flux [Jy]")
            ax3.set_title("Instrumental Pol Fluxes vs Frequency")
            leg = ax3.legend(
                frameon=True, markerscale=1, handlelength=1, loc="upper right"
            )
            leg.get_frame().set_facecolor("white")
            for le in leg.legendHandles:
                le.set_alpha(1)

            # Plot FDF
            ax4 = plt.subplot(224)
            ax4.set_xlim([-50, 50])
            ax4.plot(phi, np.abs(fdf), label=r"FDF", zorder=3)
            ax4.set_title(r"FDF : $\phi$=20 rad m$^{-2}$")
            ax4.set_xlabel("Faraday Depth [rad/m$^2$]")
            ax4.set_ylabel("Polarized Flux Density [Jy/PSF/RMSF]")
            leg = ax4.legend(
                frameon=True, markerscale=1, handlelength=1, loc="upper right"
            )
            leg.get_frame().set_facecolor("white")
            for le in leg.legendHandles:
                le.set_alpha(1)

            fig.suptitle(
                f"Amplitude Errors to Leakage/Mixing Term of X Dipoles : [{a:.5f}]",
                fontsize=16,
            )
            plt.savefig(
                f"../data/rm_theory/leak_amp_2/amp_{fr:04d}_{a:.5f}.png",
                bbox_inches="tight",
                dpi=200,
            )
            plt.close()

    #####################################################################
    #                                                                   #
    #          Save frames of leakage phase test animation              #
    #                                                                   #
    #####################################################################

    Save = True

    if Save:

        # Get stokes parameters as a function of frequency
        I_0, Q_0, U_0, V_0 = get_IQUV_complex(
            freqs, rm, ref_I_Jy, ref_V_Jy, SI, frac_pol, ref_chi=0.0, ref_freq=200e6
        )

        for fr, a in enumerate(np.linspace(0, np.pi/4, 512)):

            deg = np.rad2deg(a)
            print(f"{fr:04d}: {deg:7.3f}")

            G = 0.1 * np.exp(1j * a)

            # Apply some Hamaker leakage
            G_Ax = 1 + 0j
            G_Ay = 1 + 0j
            G_Bx = 1 + 0j
            G_By = 1 + 0j

            l_Ax = G
            l_Ay = 0 + 0j
            l_Bx = G
            l_By = 0 + 0j

            I, Q, U, V = rm_leakage(
                I_0, Q_0, U_0, V_0, G_Ax, G_Ay, G_Bx, G_By, l_Ax, l_Ay, l_Bx, l_By
            )

            # Convert stokes to instrumental pols
            XX, XY, YX, YY = stokes_instru(I, Q, U, V)

            # Determine FDF, RMSF
            fdf, rmsf, phi = rm_synth(freqs, Q, U, phi_lim=200, dphi=0.1)

            # Plot stokes vectors, FDF, RMSF
            plt.style.use("seaborn")
            fig = plt.figure(figsize=(14, 8))

            # Plot stokes vectors
            ax1 = plt.subplot(221)
            colors = plt.cm.Spectral([0.01, 0.14, 0.86, 0.99])
            ax1.set_ylabel("Flux [Jy]")
            ax1.set_title("Stoke Fluxes vs Frequency")
            stokes = [I, Q, U, V]
            for i, st in enumerate(["I", "Q", "U", "V"]):
                ax1.plot(freqs / 1e6, np.real(stokes[i]), color=colors[i], label=st)

            leg = ax1.legend(
                frameon=True, markerscale=1, handlelength=1, loc="upper right"
            )
            leg.get_frame().set_facecolor("white")
            for le in leg.legendHandles:
                le.set_alpha(1)

            # Plot RMSF
            ax2 = plt.subplot(222)
            ax2.set_xlim([-20, 20])
            ax2.plot(phi, np.abs(rmsf), label=r"$ \vert R \vert $", zorder=3)
            ax2.plot(phi, np.real(rmsf), label=r"$ real(R) $")
            ax2.plot(phi, np.imag(rmsf), label=r"$ imag(R) $")
            ax2.set_title("RMSF [-20, 20]")
            ax2.set_ylabel("RMSF")

            leg = ax2.legend(
                frameon=True, markerscale=1, handlelength=1, loc="upper right"
            )
            leg.get_frame().set_facecolor("white")
            for le in leg.legendHandles:
                le.set_alpha(1)

            # Plot Instrumental Polarizations
            ax3 = plt.subplot(223)
            colors = plt.cm.Spectral([0.01, 0.14, 0.86, 0.99])
            instru = [XX, YY, XY, YX]
            for i, inst in enumerate(["XX", "YY", "XY", "YX"]):
                ax3.plot(freqs / 1e6, np.real(instru[i]), color=colors[i], label=inst)

            ax3.set_xlabel("Frequency [MHz]")
            ax3.set_ylabel("Flux [Jy]")
            ax3.set_title("Instrumental Pol Fluxes vs Frequency")
            leg = ax3.legend(
                frameon=True, markerscale=1, handlelength=1, loc="upper right"
            )
            leg.get_frame().set_facecolor("white")
            for le in leg.legendHandles:
                le.set_alpha(1)

            # Plot FDF
            ax4 = plt.subplot(224)
            ax4.set_xlim([-50, 50])
            ax4.plot(phi, np.abs(fdf), label=r"FDF", zorder=3)
            ax4.set_title(r"FDF : $\phi$=20 rad m$^{-2}$")
            ax4.set_xlabel("Faraday Depth [rad/m$^2$]")
            ax4.set_ylabel("Polarized Flux Density [Jy/PSF/RMSF]")
            leg = ax4.legend(
                frameon=True, markerscale=1, handlelength=1, loc="upper right"
            )
            leg.get_frame().set_facecolor("white")
            for le in leg.legendHandles:
                le.set_alpha(1)

            fig.suptitle(
                f"Phase Errors to Leakage/Mixing Term of X Dipoles with 0.1 Amplitude: [{deg:7.3f}$^\circ$]",
                fontsize=16,
            )
            plt.savefig(
                f"../data/rm_theory/leak_phase_2/phase_{fr:04d}_{deg:7.3f}.png",
                dpi=200,
                bbox_inches="tight",
            )
            plt.close()

    #####################################################################
    #                                                                   #
    #             Multilayer Phase error plot for poster                #
    #                                                                   #
    #####################################################################

    Save = False

    if Save:

        # Get stokes parameters as a function of frequency
        I, Q, U, V = get_IQUV_complex(
            freqs, rm, ref_I_Jy, ref_V_Jy, SI, frac_pol, ref_chi=0.0, ref_freq=200e6
        )

        fdfs = []
        #  for fr, a in enumerate(np.linspace(0, np.pi / 20, 512)):
        for fr, a in enumerate(np.linspace(0, 5, 201)):

            #  deg = np.rad2deg(a)
            #  print(f"{fr:04d}: {deg:7.3f}")
            print(f"{fr:04d}: {a:7.3f}")

            G = 1.0 * np.exp(1j * np.deg2rad(a))

            # Apply some Hamaker leakage
            G_Ax = G
            G_Ay = 1 + 0j
            G_Bx = G
            G_By = 1 + 0j

            l_Ax = 0 + 0j
            l_Ay = 0 + 0j
            l_Bx = 0 + 0j
            l_By = 0 + 0j

            # I, Q, U, V = rm_leakage(
            #     I, Q, U, V, G_Ax, G_Ay, G_Bx, G_By, l_Ax, l_Ay, l_Bx, l_By
            # )
            I, Q, U, V = rm_leakage(
                I, Q, U, V, G_Ax, G_Ay, G_Bx, G_By, l_Ax, l_Ay, l_Bx, l_By
            )

            # Convert stokes to instrumental pols
            # XX, XY, YX, YY = stokes_instru(I, Q, U, V)

            # Determine FDF, RMSF
            fdf, rmsf, phi = rm_synth(freqs, Q, U, phi_lim=200, dphi=0.1)

            if a in np.linspace(0, 3, 4):

                fdfs.append(fdf)

            # Plot stokes vectors, FDF, RMSF
        #  plt.style.use("seaborn")

        plt.rcParams.update(
            {
                #  "font.size": 24,
                "figure.titlesize": 24,
                "axes.labelsize": 16,
                "axes.titlesize": 16,
                "xtick.labelsize": 14,
                "ytick.labelsize": 14,
                "legend.fontsize": 14,
                "font.family": "serif",
                "font.serif": "Times New Roman",
                "text.usetex": True,
            }
        )

        fig = plt.figure(figsize=(14, 5))

        # # Plot stokes vectors
        # ax1 = plt.subplot(221)
        # colors = plt.cm.Spectral([0.01, 0.14, 0.86, 0.99])
        # ax1.set_ylabel("Flux [Jy]")
        # ax1.set_title("Stoke Fluxes vs Frequency")
        # stokes = [I, Q, U, V]
        # for i, st in enumerate(["I", "Q", "U", "V"]):
        #     ax1.plot(freqs / 1e6, np.real(stokes[i]), color=colors[i], label=st)

        # leg = ax1.legend(
        #     frameon=True, markerscale=1, handlelength=1, loc="upper right"
        # )
        # leg.get_frame().set_facecolor("white")
        # for le in leg.legendHandles:
        #     le.set_alpha(1)

        # # Plot RMSF
        # ax2 = plt.subplot(222)
        # ax2.set_xlim([-20, 20])
        # ax2.plot(phi, np.abs(rmsf), label=r"$ \vert R \vert $", zorder=3)
        # ax2.plot(phi, np.real(rmsf), label=r"$ real(R) $")
        # ax2.plot(phi, np.imag(rmsf), label=r"$ imag(R) $")
        # ax2.set_title("RMSF [-20, 20]")
        # ax2.set_ylabel("RMSF")

        # leg = ax2.legend(
        #     frameon=True, markerscale=1, handlelength=1, loc="upper right"
        # )
        # leg.get_frame().set_facecolor("white")
        # for le in leg.legendHandles:
        #     le.set_alpha(1)

        # # Plot Instrumental Polarizations
        # ax3 = plt.subplot(223)
        # colors = plt.cm.Spectral([0.01, 0.14, 0.86, 0.99])
        # instru = [XX, YY, XY, YX]
        # for i, inst in enumerate(["XX", "YY", "XY", "YX"]):
        #     ax3.plot(freqs / 1e6, np.real(instru[i]), color=colors[i], label=inst)

        # ax3.set_xlabel("Frequency [MHz]")
        # ax3.set_ylabel("Flux [Jy]")
        # ax3.set_title("Instrumental Pol Fluxes vs Frequency")
        # leg = ax3.legend(
        #     frameon=True, markerscale=1, handlelength=1, loc="upper right"
        # )
        # leg.get_frame().set_facecolor("white")
        # for le in leg.legendHandles:
        #     le.set_alpha(1)

        # Plot FDF
        # ax4 = plt.subplot(224)
        ax4 = plt.subplot(111)
        ax4.set_xlim([-50, 50])

        colors = plt.cm.Spectral([0.16, 0.33, 0.83, 1])

        ls = ["solid", "dotted", "dashed", "dashdot"]

        for i in range(len(fdfs)):
            ax4.plot(
                phi,
                np.abs(fdfs[i]),
                label=rf"$\Delta\phi$ : {np.linspace(0, 3, 4)[i]}$^\circ$",
                #  zorder=-1 * i,
                linewidth=2,
                linestyle=ls[i],
                color=colors[i],
                alpha=0.9,
            )
        ax4.set_title(r"True FDF : $\phi$=20 rad m$^{-2}$")
        ax4.set_xlabel("Faraday Depth [rad/m$^2$]")
        ax4.set_ylabel("Polarized Flux Density [Jy/PSF/RMSF]")
        leg = ax4.legend(frameon=True, markerscale=1, handlelength=1, loc="upper right")
        leg.get_frame().set_facecolor("white")
        for le in leg.legendHandles:
            le.set_alpha(1)

        plt.savefig(
            #  f"../data/rm_theory/gain_phase/phase_{fr:04d}_{deg:7.3f}.png",
            f"../data/phase_sim/phase_layers_6.png",
            dpi=400,
            bbox_inches="tight",
            transparent=True,
        )
        #  plt.close()
        #  break
