import numpy as np
from matplotlib import pyplot as plt

#  plt.style.use("seaborn")
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

fig, ax = plt.subplots(figsize=(3.6, 3.6))

cols = plt.cm.Spectral(np.linspace(0, 1, 16))
x = np.arange(0.2, 1.0, 0.2)
y = np.arange(0.2, 1.0, 0.2)[::-1]

xx, yy = np.meshgrid(x, y)

left, bottom, width, height = (0.1, 0.1, 0.8, 0.8)
rect = plt.Rectangle(
    (left, bottom),
    width,
    height,
    facecolor="#eaeaf2",
    edgecolor="#DAD0C2",
    alpha=0.8,
    zorder=-1,
)
ax.add_patch(rect)
for i in range(16):
    plt.scatter(
        xx.flatten()[i],
        yy.flatten()[i],
        color=cols[i],
        edgecolor="k",
        linewidth=0.5,
        marker="s",
        s=64,
    )
    plt.text(xx.flatten()[i] + 0.04, yy.flatten()[i] + 0.04, i, fontsize=12)

plt.text(0.48, 0.95, r"$\mathcal{N}$", fontsize=18)
plt.text(0.48, -0.05, r"$\mathcal{S}$", fontsize=18)
plt.text(0.95, 0.5, r"$\mathcal{E}$", fontsize=18)
plt.text(-0.05, 0.5, r"$\mathcal{W}$", fontsize=18)

plt.axis("off")
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.tight_layout()
plt.savefig("./plots/dipole_order.pdf", bbox_inches="tight", dpi=300)
#  plt.show()
