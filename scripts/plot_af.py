"""Plot Analytic vs FEE RM Spectra."""

from matplotlib import pyplot as plt

from plot_rm import pogs_obj_loc, read_noise, read_rm_cube


if __name__=="__main__":

    pogs_fits = "../data/catalogs/POGS-II_ExGal.fits"

    pogs_obsids = {
        "1120300232" : ["POGSII-EG-337", "POGSII-EG-338", "POGSII-EG-339", "POGSII-EG-321", "POGSII-EG-313"],
        "1120300352" : ["POGSII-EG-337", "POGSII-EG-338", "POGSII-EG-339", "POGSII-EG-321", "POGSII-EG-313"],
        "1120082744" : ["POGSII-EG-034", "POGSII-EG-032", "POGSII-EG-048", "POGSII-EG-009", "POGSII-EG-010"],
        "1120082864" : ["POGSII-EG-034", "POGSII-EG-032", "POGSII-EG-048", "POGSII-EG-009", "POGSII-EG-010"],
    }

    obsids = list(pogs_obsids.keys())

    low_band = ["1120300232", "1120082744"]
    high_band = ["1120300352", "1120082864"]

    colors = ["#cc561e", "#008891"]

    for o in obsids:

        pogs = pogs_obsids[o]

        for p in pogs:

            print(f" ** INFO: Plotting RM Spec {o}: {p}")

            if o in low_band:
                title = f"{o}: ANA vs FEE: 167-200 MHz: {p}"
                fname = f"{o}_ana_fee_{p}_167-200MHz_rm.png"
            else:
                title = f"{o}: ANA vs FEE: 200-230 MHz: {p}"
                fname = f"{o}_ana_fee_{p}_200-230MHz_rm.png"

            pogs_pos = pogs_obj_loc(p, pogs_fits)

            ana_pre = f"ana_wide_{o}_"
            ana_dir = f"./data/{o}/ana_wide/imgs/cubes/"

            ra_a, dec_a, phi_a, phi_a, data_a, wcs_a = read_rm_cube(pogs_pos, f"ana_wide_{o}_", f"../data/{o}/ana_wide/imgs/cubes")
            n_a = read_noise(pogs_pos, f"ana_wide_{o}_", f"../data/{o}/ana_wide/imgs/cubes")

            ra_f, dec_f, phi_f, phi_f, data_f, wcs_f = read_rm_cube(pogs_pos, f"fee_wide_{o}_", f"../data/{o}/fee_wide/imgs/cubes")
            n_f = read_noise(pogs_pos, f"fee_wide_{o}_", f"../data/{o}/fee_wide/imgs/cubes")


            # Plotting stuff
            plt.style.use("seaborn")

            fig = plt.figure(figsize=(10, 7))
            ax = fig.add_subplot(1, 1, 1)
            ax.plot(
                phi_f,
                data_f[dec_f, ra_f, :] - n_f,
                label="FEE",
                linewidth=2,
                color=colors[0],
                alpha=0.9,
            )
            ax.plot(
                phi_a,
                data_a[dec_a, ra_a, :] - n_a,
                label="ANA",
                linewidth=2,
                color=colors[1],
                alpha=0.9,
            )

            ax.set_xlabel("Faraday Depth [rad/m$^2$]")
            ax.set_ylabel("Polarized Flux Density [Jy/PSF/RMSF]")
            ax.set_title(title)
            ax.grid(True, color="white", linewidth=1.2, alpha=0.9, ls="dotted")

            leg = plt.legend(frameon=True, markerscale=1, handlelength=1)
            leg.get_frame().set_edgecolor("#22222233")
            leg.get_frame().set_facecolor("none")
            for le in leg.legendHandles:
                le.set_alpha(1)

            plt.tight_layout()
            plt.savefig(f"../data/leakage_plots/{fname}", dpi=300)
