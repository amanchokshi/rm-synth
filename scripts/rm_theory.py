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


def leakage(XX, XY, YX, YY, D):
    """Hamaker Dipole leakage."""

    XX_measured = XX + D * XY + D * YX + D ** 2 * YY
    XY_measured = -D * XX + XY - D ** 2 * YX + D * YY
    YX_measured = -D * XX - D ** 2 * XY + YX + D * YY
    YY_measured = D ** 2 * XX - D * XY - D * YX + YY

    return XX_measured, XY_measured, YX_measured, YY_measured


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

    # Convert stokes to instrumental pols
    XX, XY, YX, YY = stokes_instru(I, Q, U, V)
    print(XX[0])

    # Leakage as defined by Hamaker - I
    XX, XY, YX, YY = leakage(XX, XY, YX, YY, 0.0)
    print(XX[0])

    # Convert back to stokes
    I, Q, U, V = instru_stokes(XX, XY, YX, YY)

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

    leg = ax1.legend(frameon=True, markerscale=1, handlelength=1)
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

    leg = ax2.legend(frameon=True, markerscale=1, handlelength=1)
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
    leg = ax3.legend(frameon=True, markerscale=1, handlelength=1)
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    # Plot FDF
    ax4 = plt.subplot(224)
    ax4.set_xlim([-10, 50])
    ax4.plot(phi, np.abs(fdf), label=r"FDF", zorder=3)
    ax4.set_title(r"FDF : $\phi$=20 rad m$^{-2}$")
    ax4.set_xlabel("Faraday Depth [rad/m$^2$]")
    ax4.set_ylabel("Polarized Flux Density [Jy/PSF/RMSF]")
    leg = ax4.legend(frameon=True, markerscale=1, handlelength=1)
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    plt.tight_layout()
    plt.show()
