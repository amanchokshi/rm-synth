import corner
import emcee
import numpy as np
from chainconsumer import ChainConsumer
from matplotlib import pyplot as plt

if __name__ == "__main__":

    reader = emcee.backends.HDFBackend("./S06YY_rf1YY_825491_beam_mcmc.h5")
    #  samples = reader.get_chain(discard=5491)
    #  samples = reader.get_chain()
    #  samples = samples[7000:7100, :, :]
    #  print(samples.shape)

    # fig, axes = plt.subplots(2, figsize=(10, 7), sharex=True)
    # #  samples = sampler.get_chain()
    # labels = ["$A_1$", "$A_2$"]
    # for i in range(2):
    #     ax = axes[i]
    #     ax.plot(samples[:, :, i], "k", alpha=0.1)
    #     ax.set_xlim(0, len(samples))
    #     ax.set_ylabel(labels[i])
    #     #  ax.yaxis.set_label_coords(-0.1, 0.5)

    # axes[-1].set_xlabel("iteration")
    # #  plt.tight_layout()
    # plt.savefig("convergence_2_amps_v3.png")
    # plt.close()

    # flat_samples = reader.get_chain(flat=True)
    # # fig = corner.corner(flat_samples, show_titles=True)
    # fig = corner.corner(flat_samples, labels=["A1", "A2"], quantiles=[0.16, 0.5, 0.84], truths=[0.5, 0.8], show_titles=True)
    # plt.savefig("emcee_2_amps_v2.png")
    # plt.close()
    #  colors = spec(np.linspace(0, 1, 13))

    #  posterior = reader.get_log_prob(flat=True)

    samples = reader.get_chain(discard=5491, flat=True)
    d_3 = samples[:, 3]
    d_10 = samples[:, 10]

    ss = np.array([d_3, d_10]).T

    ndim = 2
    c = ChainConsumer()
    c.add_chain(
        ss,
        parameters=[
            #  "$d_0$",
            #  "$d_1$",
            #  "$d_2$",
            "$d_3$",
            #  "$d_4$",
            #  "$d_5$",
            #  "$d_6$",
            #  "$d_7$",
            #  "$d_8$",
            #  "$d_9$",
            "$d_{10}$",
            #  "$d_{11}$",
            #  "$d_{12}$",
            #  "$d_{13}$",
            #  "$d_{14}$",
            #  "$d_{15}$",
        ],
        name="samples",
        walkers=1,
    )
    c.configure(
        statistics="cumulative",
        max_ticks=10,
        sigma2d=True,
        diagonal_tick_labels=False,
        colors=["#5E4FA1"],
        plot_hists=False,
        tick_font_size=6,
        spacing=1.4,
        linewidths=0.4,
    )
    c.configure_truth(
        linewidth=1.4, linestyle="solid", color="#FCAD61",
    )
    S06YY_T = np.array(
        [
            #  0.8671226184660479,
            #  0.7816316560820713,
            #  0.9956033219345384,
            0.8466047874938935,
            #  1.0,
            #  0.885686370297997,
            #  0.9863214460185638,
            #  0.6893014167073767,
            #  1.0,
            #  0.8382999511480215,
            1.0,
            #  0.761113825109917,
            #  0.986809965803615,
            #  0.6082071323888617,
            #  1.0,
            #  0.9848558866634098,
        ]
    )
    plt.style.use("seaborn")
    plt.rcParams.update(
        {
            "text.usetex": True,
            "font.family": "serif",
            "axes.labelsize": 7,
            "axes.titlesize": 9,
            "font.size": 8,
            "legend.fontsize": 10,
        }
    )

    #  fig = plt.figure(figsize=(7.6, 9.4))
    #  fig = c.plotter.plot(figsize=(20, 20), truth=S06YY_T)
    fig = c.plotter.plot(figsize=(1.7, 1.7), truth=S06YY_T)
    plt.savefig("S06YY_emcee_82000_truths_3_10.pdf", bbox_inches="tight", dpi=600)
