from plot_rm import pogs_obj_loc, read_rm_cube

if __name__ == "__main__":

    from matplotlib import pyplot as plt
    import matplotlib.pylab as pl
    import numpy as np

    pogs_pos = pogs_obj_loc("POGSII-EG-321", "../slurm/iono/POGS-II_ExGal.fits")

    patches = ["fee_pa10", "fee_pa33", "fee_pa66", "fee_pa100", "fee_pa333", "fee_pa666", "fee_pa1000"]

    plt.style.use("seaborn")

    colors = pl.cm.Spectral_r(np.linspace(0,1,len(patches)))

    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(1, 1, 1)

    for i, patch in enumerate(patches):
        ra_x, dec_y, phi_z, phi, data_p, wcs = read_rm_cube(pogs_pos, "rts_imgr_", f"../data/1086351512/{patch}/run_i/cubes_iono")
        ax.plot(phi, data_p[dec_y, ra_x, :], label=f"{patch}", color=colors[i], alpha=0.9)

    ax.set_xlabel("Faraday Depth [rad/m$^2$]")
    ax.set_ylabel("Polarized Flux Density [Jy/PSF/RMSF]")
    ax.set_title("POGSII-EG-321 RTS Patch RM effects")

    leg = plt.legend(frameon=True, markerscale=1, handlelength=1)
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    plt.tight_layout()
    plt.savefig(f"rm_patch_test.png", dpi=300)
