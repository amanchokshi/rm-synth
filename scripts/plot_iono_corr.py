from plot_rm import pogs_obj_loc, read_rm_cube

if __name__ == "__main__":

    from matplotlib import pyplot as plt

    #####################################################################
    #                                                                   #
    #                   Read in data and create spectra                 #
    #                                                                   #
    #####################################################################

    pogs_pos = pogs_obj_loc("POGSII-EG-321", "../slurm/iono/POGS-II_ExGal.fits")

    ra_x, dec_y, phi_z, phi, data_p, wcs = read_rm_cube(pogs_pos, "rts_imgr_", "../data/1086351512/rts_imgr/run_i/cubes")
    ra_x_i, dec_y_i, phi_z_i, phi_i, data_p_i, wcs_i = read_rm_cube(pogs_pos, "rts_imgr_", "../data/1086351512/rts_imgr/run_i/cubes_iono")

    # Elegant Seaborn
    plt.style.use("seaborn")

    #####################################################################
    #                                                                   #
    #                  Plot RM spectra of POGS source                   #
    #                                                                   #
    #####################################################################

    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(phi, data_p[dec_y, ra_x, :], label="Raw RM", color="#822659", alpha=0.8)
    ax.plot(phi_i, data_p_i[dec_y_i, ra_x_i, :], label="Iono Corr RM", color="#207561")
    ax.set_xlabel("Faraday Depth [rad/m$^2$]")
    ax.set_ylabel("Polarized Flux Density [Jy/PSF/RMSF]")
    ax.set_title("POGSII-EG-321 Ionospheric RM effects")

    leg = plt.legend(frameon=True, markerscale=1, handlelength=1)
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    plt.tight_layout()
    plt.savefig(f"rm_spec.png", dpi=300)
