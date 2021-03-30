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


def rm_synth(freqs, Q, U, phi_lim=200, dphi=0.5):
    """Do RM Synthesis on stokes Q & U vectors."""

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
    I, Q, U, V = get_IQUV_complex(
        freqs, rm, ref_I_Jy, ref_V_Jy, SI, frac_pol, ref_chi=0.0, ref_freq=200e6
    )

    #  Q = 0.1 * I + Q
    #  I = 0.9 * I

    # Determine FDF, RMSF
    fdf, rmsf, phi = rm_synth(freqs, Q, U, phi_lim=200, dphi=0.1)

    # Plot stokes vectors, FDF, RMSF
    plt.style.use("seaborn")
    fig = plt.figure(figsize=(10, 7))

    ax1 = plt.subplot(121)
    colors = plt.cm.Spectral([0.01, 0.14, 0.86, 0.99])

    stokes = [I, Q, U, V]

    for i, st in enumerate(["I", "Q", "U", "V"]):
        ax1.plot(freqs / 1e6, np.real(stokes[i]), color=colors[i], label=st)

    ax1.set_xlabel("Frequency [MHz]")
    ax1.set_ylabel("Stokes Flux [Jy]")
    ax1.set_title("Stoke Fluxes vs Frequency")

    leg = ax1.legend(frameon=True, markerscale=1, handlelength=1)
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    ax2 = plt.subplot(222)
    ax2.plot(phi, np.abs(rmsf), label=r"$ \vert R \vert $", zorder=3)
    ax2.plot(phi, np.real(rmsf), label=r"$ real(R) $")
    ax2.plot(phi, np.imag(rmsf), label=r"$ imag(R) $")
    ax2.set_xlim([-20, 20])

    leg = ax2.legend(frameon=True, markerscale=1, handlelength=1)
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    ax2.set_title("RMSF [-20, 20]")
    #  ax2.set_xlabel("Faraday Depth [rad/m$^2$]")
    ax2.set_ylabel("RMSF")

    ax3 = plt.subplot(224)
    ax3.plot(phi, np.abs(fdf), label=r"FDF", zorder=3)
    ax3.set_xlim([-10, 50])

    ax3.set_title(r"FDF : $\phi$=20 rad m$^{-2}$")
    ax3.set_xlabel("Faraday Depth [rad/m$^2$]")
    ax3.set_ylabel("Polarized Flux Density [Jy/PSF/RMSF]")

    leg = ax3.legend(frameon=True, markerscale=1, handlelength=1)
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    plt.tight_layout()
    plt.show()
