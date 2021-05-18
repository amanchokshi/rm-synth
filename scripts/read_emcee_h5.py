import corner
import emcee
import numpy as np
from matplotlib import pyplot as plt

if __name__ == "__main__":

    reader = emcee.backends.HDFBackend("./beam_mcmc.h5")
    samples = reader.get_chain()

    # fig, axes = plt.subplots(16, figsize=(10, 7), sharex=True)
    # #  samples = sampler.get_chain()
    # #  labels = ["A", r"$\beta$", "B", r"$\omega$"]
    # for i in range(16):
    #     ax = axes[i]
    #     ax.plot(samples[:, :, i], "k", alpha=0.2)
    #     ax.set_xlim(0, len(samples))
    #     #  ax.set_ylabel(labels[i])
    #     #  ax.yaxis.set_label_coords(-0.1, 0.5)

    # axes[-1].set_xlabel("iteration")
    # plt.show()

    flat_samples = reader.get_chain(flat=True)
    fig = corner.corner(flat_samples, quantiles=[0.16, 0.5, 0.84], show_titles=True)
    plt.savefig("emcee.png")
