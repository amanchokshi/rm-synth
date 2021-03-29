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


def get_IQU_complex(freqs, rm, ref_I_Jy, SI, frac_pol, ref_chi=0.0, ref_freq=200e6):
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

    return I, Q, U


def plot_iquv(freqs, I, Q, U, V):
    """Plot Stokes parameters vs frequency.

    Parameters
    ----------
    freqs : numpy.array / float
        An array of freqencies in Hz at which flux is to be determined
    I, Q, U, V : numpy.array
        Arrays of complex stokes fluxes

    Returns
    -------
    matplotlib.pyplot.figure
    """

    plt.style.use("seaborn")
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111)

    colors = plt.cm.Spectral([0.01, 0.14, 0.86, 0.99])

    stokes = [I, Q, U, V]

    for i, st in enumerate(["I", "Q", "U", "V"]):
        plt.plot(freqs / 1e6, np.real(stokes[i]), color=colors[i], label=st)

    ax.set_xlabel("Frequency [MHz]")
    ax.set_ylabel("Stokes Flux [Jy]")
    ax.set_title("Stoke Fluxes vs Frequency")

    leg = plt.legend(frameon=True, markerscale=1, handlelength=1)
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
    rm = 49
    frac_pol = 0.21
    ref_I_Jy = 7
    ref_V_Jy = 1

    # Arrays of freqencies in Hz
    freqs = np.arange(low_freq, high_freq + fine_channel, fine_channel)

    # Get stokes parameters as a function of frequency
    I, Q, U = get_IQU_complex(
        freqs, rm, ref_I_Jy, SI, frac_pol, ref_chi=0.0, ref_freq=200e6
    )
    V = spectal_index(freqs, SI, ref_V_Jy, ref_freq=200e6)

    # Complex linear polarization
    P = Q + 1j * U

    # Plot I, Q, U, V
    plot_iquv(freqs, I, Q, U, V)
