"""
Create a RTS style source list for a
polarized source with given
rotation measure
"""

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib import pyplot as plt

from rm_theory import get_IQUV_complex, rm_synth

if __name__ == "__main__":

    # MWA constants
    low_freq = 168e6
    high_freq = 232e6
    fine_channel = 5e3

    # RM Source Constants
    SI = -0.7
    rm = 20
    frac_pol = 0.20
    ref_I_Jy = 10
    ref_V_Jy = 1

    # Arrays of freqencies in Hz
    freqs = np.arange(low_freq, high_freq + fine_channel, fine_channel)

    # Get stokes parameters as a function of frequency
    I, Q, U, V = get_IQUV_complex(
        freqs, rm, ref_I_Jy, ref_V_Jy, SI, frac_pol, ref_chi=0.0, ref_freq=200e6
    )

    rm_coords = SkyCoord(
        ra=10.44137913863797 * u.degree, dec=-26.78792973842179 * u.degree, frame="icrs"
    )

    # with open("../data/rts_srclists/srclist_eor1_rm001_v2.txt", "w") as f:
    #     f.write(f"SOURCE RM-001 {rm_coords.ra.hour} {rm_coords.dec.deg}\n")

    #     for i, freq in enumerate(freqs):
    #         f.write(
    #             f"FREQ {freq:.6e} {I[i].real:.5f} {Q[i].real:.5f} {U[i].real:.5f} {V[i].real:.5f}\n"
    #         )
    #     f.write("ENDSOURCE\n")

    FDF, RMSF, phi_arr = rm_synth(freqs, Q, U, phi_lim=200, dphi=0.1)

    plt.style.use("seaborn")
    plt.plot(phi_arr, np.abs(FDF))
    plt.tight_layout()
    plt.show()
    # plt.plot(freqs, I.real, label="I")
    # plt.plot(freqs, Q.real, label="Q")
    # plt.plot(freqs, U.real, label="U")
    # plt.plot(freqs, V.real, label="V")
    # leg = plt.legend(loc="upper right", frameon=True)
    # leg.get_frame().set_facecolor("white")
    # for le in leg.legendHandles:
    #     le.set_alpha(1)
    # plt.ylabel("FLux [Jy]")
    # plt.xlabel("Frequency [Hz]")

    # plt.plot(phi_arr, np.abs(FDF))
    # plt.xlim([-50, 50])

    # plt.tight_layout()
    # plt.show()
