"""Monte Carle to determine MWA beam dipole amplitudes."""

import mwa_hyperbeam
import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm

import healpix as hpx


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


def chisq_prob(data=None, model=None):
    """Chi-square goodness of fit test.

    Determine the probability that a model fits the data

    Parameter
    ---------
    data : numpy.array
        A data array
    model : numpy.array
        A model array of same size as data

    Returns
    -------
    p_value : float
        The probability that the model fits the data
    """
    chisq = np.sum(np.square(data - model) / model)
    p_value = np.exp(-0.5 * chisq)

    return p_value, chisq


if __name__ == "__main__":

    # We can make a new beam object with a path to the HDF5 file specified by MWA_BEAM_FILE.
    beam = mwa_hyperbeam.FEEBeam()

    # Healpix map with given nside
    nside = 32

    # Zenith angle and Azimuth of healpix pixels
    za, az = hpx.healpix_za_az(nside=nside)

    # Satellite beam map frequency
    freq = 138e6

    # Zenith Pointing
    delays = [0] * 16

    # Normalize maps to zenith
    norm_to_zenith = True

    # Setup a beam with the first dipole dead
    # We'll try and determine which dipole is dead using MCMC magic
    # amps_15 = [0.0] * 1 + [1.0] * 15
    # jones_15 = beam.calc_jones_array(az, za, freq, delays, amps_15, norm_to_zenith)
    # unpol_beam_15 = makeUnpolInstrumentalResponse(jones_15, jones_15)

    # XX_15 = np.real(unpol_beam_15[:, 0])
    # YY_15 = np.real(unpol_beam_15[:, 3])

    # Monte Carlo Magic Below
    # define a set of N random dipole amps between 0, 1

    # N = 1000

    # # Array of uniformly distributed amplitudes
    # # The first column represents the initial conditions of the MCMC
    # amps_rand = np.random.uniform(low=0.0, high=1.0, size=(16, N + 1))

    # # Initialisation
    # amps_ini = amps_rand[:, 0]

    # # Evaluate the beam with the initial amplitudes
    # jones = beam.calc_jones_array(az, za, freq, delays, amps_ini, norm_to_zenith)
    # unpol_beam = makeUnpolInstrumentalResponse(jones, jones)

    # XX = np.real(unpol_beam[:, 0])
    # #  YY = np.real(unpol_beam[:, 3])

    # # This is the probability of the initial amplitude guess
    # prob_old, _ = chisq_prob(data=XX_15, model=XX)

    # amps_accepted = [amps_ini]
    # chi_sq = []

    # for i in tqdm(range(N)):

    #     # A candidate jump
    #     amps = amps_rand[:, i + 1]

    #     jones = beam.calc_jones_array(az, za, freq, delays, amps, norm_to_zenith)
    #     unpol_beam = makeUnpolInstrumentalResponse(jones, jones)

    #     XX = np.real(unpol_beam[:, 0])
    #     #  YY = np.real(unpol_beam[:, 3])

    #     # Probability of new model, with random amps
    #     prob_new, chisq = chisq_prob(data=XX_15, model=XX)

    #     # Evaluate whether it's worth making this jump
    #     if prob_new / prob_old > np.random.rand():

    #         # The new probability becomes the initial condition for the next loop
    #         prob_old = prob_new
    #         amps_accepted.append(amps)
    #         chi_sq.append(chisq)

    # amps_accepted = np.array(amps_accepted)
    # chi_sq = np.array(chi_sq)

    # np.save("./mcmc/amps_accepted.npy", amps_accepted)
    # np.save("./mcmc/chisq_accepted.npy", chi_sq)

    # MCMC Below

    # Number of jumps to try
    N = 300000

    # Standard deviation of normal distribution from which to draw candidate jumps
    #  step = 0.14
    step = 0.2

    # Lets create synthetic data at try to recover the input parameters
    amps_15 = [0.0] * 1 + [1.0] * 15
    jones_15 = beam.calc_jones_array(az, za, freq, delays, amps_15, norm_to_zenith)
    unpol_beam_15 = makeUnpolInstrumentalResponse(jones_15, jones_15)

    XX_15 = np.real(unpol_beam_15[:, 0])
    YY_15 = np.real(unpol_beam_15[:, 3])

    # Limits of dipole amplitudes
    amp_min = 0
    amp_max = 1

    amps_proposal = np.random.uniform(low=amp_min, high=amp_max, size=16)

    # Evaluate the beam with the initial amplitudes
    jones = beam.calc_jones_array(az, za, freq, delays, amps_proposal, norm_to_zenith)
    unpol_beam = makeUnpolInstrumentalResponse(jones, jones)

    XX = np.real(unpol_beam[:, 0])
    #  YY = np.real(unpol_beam[:, 3])

    # This is the probability of the initial amplitude guess
    prob_old, chi_old = chisq_prob(data=XX_15, model=XX)

    amps_old = amps_proposal

    amps_list = [amps_old]
    chisq_list = [chi_old]

    # Jump proposal

    count = 0

    while count <= N:

        #  print(count)

        jump = False
        while jump is False:

            amps_proposal = np.array(
                [
                    np.random.normal(loc=amps_old[i], scale=step, size=1)[0]
                    for i in range(16)
                ]
            )

            if np.all(amps_proposal <= amp_max) & np.all(amps_proposal >= amp_min):
                jump = True
                count += 1

                # Evaluate the beam with proposed amplitudes
                jones = beam.calc_jones_array(
                    az, za, freq, delays, amps_proposal, norm_to_zenith
                )
                unpol_beam = makeUnpolInstrumentalResponse(jones, jones)

                XX = np.real(unpol_beam[:, 0])
                #  YY = np.real(unpol_beam[:, 3])

                # This is the probability of the initial amplitude guess
                prob_new, chi_new = chisq_prob(data=XX_15, model=XX)

                if (prob_new / prob_old) > np.random.rand():
                    amps_list.append(amps_proposal)
                    chisq_list.append(chi_new)

                    # This successful jump becomes the initial condition
                    # of the next jump - a markov chain has been formed
                    prob_old = prob_new
                    amps_new = amps_proposal

    amps_list = np.array(amps_list)
    chisq_list = np.array(chisq_list)

    np.save("../data/mcmc/amps_list.npy", amps_list)
    np.save("../data/mcmc/chisq_list.npy", chisq_list)

    #  plt.style.use("seaborn")
    #  plt.plot(range(len(chisq_list)), chisq_list, "-o", color="seagreen")
    #  plt.xlabel("Jump Number")
    #  plt.ylabel("Chi Square")
    #  plt.tight_layout()
    #  plt.show()

    #  plt.hist(np.array(amps_list)[:, 0], bins=30, ec="k")
    #  plt.tight_layout()
    #  plt.show()
