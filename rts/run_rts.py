#!/usr/bin/env python

def make_rts_setup(
    obs=None,
    out_dir=None,
    tag=None,
    data_dir=None,
    gpu_box=None,
    fscrunch=None,
    phase_ra=None,
    phase_dec=None,
):
    """Create the rts_setup script for the given observation."""

    # Default source list
    srclist_dir = "/pawsey/mwa/software/python3/srclists/master"
    srclist_file = "srclist_pumav3_EoR0aegean_EoR1pietro+ForA"

    with open(f"{out_dir}/{obs}/{tag}/rts_setup_{obs}.sh", "w+") as outfile:

        outfile.write("#!/bin/bash -l\n")
        outfile.write("\n")

        outfile.write(f'#SBATCH --job-name="set_{obs}"\n')
        outfile.write(f"#SBATCH --output=RTS-setup-{obs}-%A.out\n")
        outfile.write(f"#SBATCH --error=RTS-setup-{obs}-%A.err\n")
        outfile.write("#SBATCH --nodes=1\n")
        outfile.write("#SBATCH --ntasks-per-node=1\n")
        outfile.write("#SBATCH --time=00:20:00\n")
        outfile.write("#SBATCH --clusters=garrawarla\n")
        outfile.write("#SBATCH --partition=gpuq\n")
        outfile.write("#SBATCH --account=mwaeor\n")
        outfile.write("#SBATCH --gres=gpu:1\n")
        outfile.write("\n")

        outfile.write("module use /pawsey/mwa/software/python3/modulefiles\n")
        outfile.write("module load python-singularity\n")
        outfile.write("module load srclists/master\n")
        outfile.write("module load mongoose\n")
        outfile.write("module load mwa_pb\n")
        outfile.write("\n")

        outfile.write("set -eux\n")
        outfile.write("\n")

        outfile.write("# Create symlinks from gpubox files to out_dir\n")
        outfile.write(f"cd {out_dir}/{obs}/{tag}\n")
        outfile.write(f"ln -s {data_dir}/{obs}/{gpu_box}/* .\n")
        outfile.write("\n")

        outfile.write("# Generate RTS mwaf files.\n")
        outfile.write("reflag-mwaf-files\n")
        outfile.write("\n")

        outfile.write("# Generate a source list for the patch step.\n")
        outfile.write("srclist_by_beam.py -n 1000 \\\n")
        outfile.write(
            f"                   --srclist={srclist_dir}/{srclist_file}.txt \\\n"
        )
        outfile.write(
            f"                   --metafits={out_dir}/{obs}/{tag}/{obs}.metafits\n"
        )
        outfile.write("\n")

        outfile.write("# Generate a source list for the peel step.\n")
        outfile.write("srclist_by_beam.py -n 3000 \\\n")
        outfile.write(
            f"                   --srclist={srclist_dir}/{srclist_file}.txt \\\n"
        )
        outfile.write(
            f"                   --metafits={out_dir}/{obs}/{tag}/{obs}.metafits \\\n"
        )
        outfile.write("                   --no_patch \\\n")
        outfile.write("                   --cutoff=90\n")
        outfile.write("\n")

        outfile.write("# Generate the RTS .in files for both patching and peeling.\n")
        outfile.write("rts-in-file-generator patch \\\n")
        outfile.write(f"                      --base-dir {out_dir}/{obs}/{tag} \\\n")

        if phase_ra and phase_dec:
            outfile.write(f"                      --force-ra {phase_ra} \\\n")
            outfile.write(f"                      --force-dec {phase_dec} \\\n")

        outfile.write(f"                      --fscrunch {fscrunch} \\\n")
        outfile.write(
            f"                      --metafits {out_dir}/{obs}/{tag}/{obs}.metafits \\\n"
        )
        outfile.write(
            f"                      --srclist {srclist_file}_{obs}_patch1000.txt \\\n"
        )
        outfile.write("                      --num-primary-cals 1 \\\n")
        outfile.write("                      > rts_patch.in\n")
        outfile.write("\n")

        outfile.write("rts-in-file-generator peel \\\n")
        outfile.write(f"                      --base-dir {out_dir}/{obs}/{tag} \\\n")

        if phase_ra and phase_dec:
            outfile.write(f"                      --force-ra {phase_ra} \\\n")
            outfile.write(f"                      --force-dec {phase_dec} \\\n")

        outfile.write(f"                      --fscrunch {fscrunch} \\\n")
        outfile.write(
            f"                      --metafits {out_dir}/{obs}/{tag}/{obs}.metafits \\\n"
        )
        outfile.write(
            f"                      --srclist {srclist_file}_{obs}_peel3000.txt \\\n"
        )
        outfile.write("                      --num-primary-cals 5 \\\n")
        outfile.write("                      --num-cals 1000 \\\n")
        outfile.write("                      --num-peel 1000 \\\n")
        outfile.write("                      > rts_peel.in\n")
        outfile.write("\n")

        outfile.write("# Ensure permissions are sensible.\n")
        outfile.write("find . -user $USER -type d -exec chmod g+rwx,o+rx,o-w {} \;\n")
        outfile.write("find . -user $USER -type f -exec chmod g+rw,o+r,o-w {} \;\n")
        outfile.write("\n")

        outfile.write("echo gator rts_setup.sh finished successfully.\n")

    return f"{out_dir}/{obs}/{tag}/rts_setup_{obs}.sh"


def make_rts_run(obs=None, out_dir=None, tag=None, clean=None):
    """Create the rts_setup script for the given observation."""

    with open(f"{out_dir}/{obs}/{tag}/rts_run_{obs}.sh", "w+") as outfile:

        outfile.write("#!/bin/bash -l\n")
        outfile.write("\n")

        outfile.write(f'#SBATCH --job-name="rts_{obs}"\n')
        outfile.write(f"#SBATCH --output=RTS-run-{obs}-%A.out\n")
        outfile.write(f"#SBATCH --error=RTS-run-{obs}-%A.err\n")
        outfile.write("#SBATCH --nodes=25\n")
        outfile.write("#SBATCH --ntasks-per-node=1\n")
        outfile.write("#SBATCH --time=00:30:00\n")
        outfile.write("#SBATCH --clusters=garrawarla\n")
        outfile.write("#SBATCH --partition=gpuq\n")
        outfile.write("#SBATCH --account=mwaeor\n")
        outfile.write("#SBATCH --gres=gpu:1\n")
        outfile.write("\n")

        outfile.write("module use /pawsey/mwa/software/python3/modulefiles\n")
        outfile.write("module load RTS/sla_to_pal\n")
        outfile.write("module load python-singularity\n")
        outfile.write("\n")

        outfile.write("set -eux\n")
        outfile.write("command -v rts_gpu\n")
        outfile.write("export UCX_MEMTYPE_CACHE=n\n")
        outfile.write("date\n")
        outfile.write("srun -n 25 --export=ALL rts_gpu rts_patch.in\n")
        outfile.write("date\n")
        outfile.write("\n")

        outfile.write("# Allow python scripts to fail.\n")
        outfile.write("set +e\n")
        outfile.write(
            "srun -n 1 -N 1 --export=ALL plot_BPcal_128T.py --both --outname BPcal.png\n"
        )
        outfile.write(
            f"srun -n 1 -N 1 --export=ALL plot_CalSols.py --base_dir=`pwd` -n {obs}\n"
        )
        outfile.write("set -e\n")
        outfile.write("date\n")
        outfile.write("\n")

        outfile.write("srun -n 25 --export=ALL rts_gpu rts_peel.in\n")
        outfile.write("date\n")
        outfile.write("\n")

        outfile.write("mkdir uvfits && mv *.uvfits uvfits\n")
        outfile.write("cp -rL *metafits*.fits uvfits\n")
        outfile.write("\n")

        outfile.write("# Ensure permissions are sensible!\n")
        outfile.write("find . -user $USER -type d -exec chmod g+rwx,o+rx,o-w {} \;\n")
        outfile.write("find . -user $USER -type f -exec chmod g+rw,o+r,o-w {} \;\n")

    return f"{out_dir}/{obs}/{tag}/rts_run_{obs}.sh"


if __name__ == "__main__":

    import argparse
    import json
    import os
    import re
    from pathlib import Path
    from subprocess import check_output
    from sys import exit

    import numpy as np

    parser = argparse.ArgumentParser(description="Setup and run RTS jobs on Garrawarla")

    parser.add_argument(
        "--tag",
        metavar="\b",
        required=True,
        help="Unique tag for each rts run, will create a <tag> subdir for rts outputs in the out_dir/obsid dir/<tag>",
    )

    parser.add_argument(
        "--rts_in",
        metavar="\b",
        default="/astro/mwaeor/achokshi/rm-synth/rts/rts_in.json",
        help="Json file containing obsids, flagged channels and tiles. Default:/astro/mwaeor/achokshi/rm-synth/rts_in.json",
    )

    parser.add_argument(
        "--out_dir",
        metavar="\b",
        default="/astro/mwaeor/achokshi/rm-synth/data",
        help="Base output dir. Default:/astro/mwaeor/achokshi/rm-synth/data",
    )

    parser.add_argument(
        "--data_dir",
        metavar="\b",
        default="/astro/mwaeor/MWA/data",
        help="Base data dir. Default:/astro/mwaeor/MWA/data",
    )

    parser.add_argument(
        "--gpu_box",
        metavar="\b",
        default="gpubox",
        help="GPU boxfiles dir within the obsid dir. Default:gpubox",
    )

    parser.add_argument(
        "--fscrunch",
        metavar="\b",
        default=1,
        type=int,
        help="The number of channels to average. Default:1",
    )

    parser.add_argument(
        "--phase_ra", metavar="\b", help="Force RA phase centre in degrees",
    )

    parser.add_argument(
        "--phase_dec", metavar="\b", help="Force DEC phase centre in degrees",
    )

    parser.add_argument(
        "--no_run",
        action="store_true",
        help="<FLAG> - Don't submit the jobs to SLURM queue after creating rts in files",
    )

    parser.add_argument(
        "--clean",
        action="store_true",
        help="<FLAG> - Clean out RTS directories of everything but metadata and uvfits, then exit",
    )

    parser.add_argument(
        "--check",
        action="store_true",
        help="<FLAG> - Run a check on previous specified run for uvfits outputs, then exit",
    )

    args = parser.parse_args()

    # Unique tag to /out_dir/obsid/<tag>
    tag = args.tag

    # Check if rts_in.json exists
    rts_in = Path(args.rts_in)
    if not rts_in.is_file():
        print("Invalid rts_in json file, enter valid path")

    # Check if out_dir exists and create if doesn't
    out_dir = Path(args.out_dir)
    if not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)
        print("Output dir did not exist. It has been created")

    # Check if data_dir exists
    data_dir = Path(args.data_dir)
    if not data_dir.exists():
        print("Data dir did not exist, enter valid path")

    # Name of dir with gpubox files
    gpu_box = args.gpu_box

    with open(rts_in) as f:
        config = json.load(f)

        try:
            obsids = config["obsids"]
        except Exception:
            print("Check rts_in json file. No obsids found ")

        # Check if flagged chans exist
        if "flagged_channels" in config.keys():
            flagged_chans = config["flagged_channels"]
        else:
            flagged_chans = None

        # Check if flagged tiles exist
        if "flagged_tiles" in config.keys():
            flagged_tiles = config["flagged_tiles"]
        else:
            flagged_tiles = None

    # If set, clean out the RTS directories and exit
    if args.clean:
        for obs in obsids:

            # Find and nuke all files in out_dir. Ignores subdirs
            os.popen(f"rm {out_dir}/{obs}/{tag}/* &> /dev/null")
            print(f"Cleaned {out_dir}/{obs}/{tag}")

    # If set, check the output of all jobs for expected files+size and exit
    elif args.check:
        for obs in obsids:

            # check that 24 uvfits exist
            uv_dir = Path(f"{out_dir}/{obs}/{tag}/uvfits/")
            num_uvfits = len([k.name for k in uv_dir.glob("*.uvfits")])

            if num_uvfits < 24:
                print(f"{obs}/{tag} does not have 24 uvfits files")
            else:
                # now check that they surpass an expected size
                size_uvfits = sum(
                    f.stat().st_size for f in uv_dir.glob("*.uvfits") if f.is_file()
                )

                # 4.2G is appropriate for 40kHz, 8 second uvfits.
                if int(size_uvfits) < 4200000/args.fscrunch:
                    print(f"{obs}/{tag} uvfits are less than 2.1G")

                else:
                    print(f"Check complete. 24 uvfits files are present in {obs}/{tag} and look good")

    else:
        # Write all shell scripts for running the rts
        setup_jobs = []
        run_jobs = []
        for obs in obsids:

            # Create rts output dir
            rts_out = Path(f"{out_dir}/{obs}/{tag}")
            rts_out.mkdir(parents=True, exist_ok=True)

            if flagged_chans:
                np.array(flagged_chans).tofile(
                    f"{rts_out}/flagged_channels.txt", "\n", "%s"
                )

            if flagged_tiles:
                np.array(flagged_tiles).tofile(
                    f"{rts_out}/flagged_tiles.txt", "\n", "%s"
                )

            setup_job = make_rts_setup(
                obs=obs,
                out_dir=out_dir,
                tag=tag,
                data_dir=data_dir,
                gpu_box=gpu_box,
                fscrunch=args.fscrunch,
                phase_ra=args.phase_ra,
                phase_dec=args.phase_dec,
            )
            setup_jobs.append(setup_job)

            run_job = make_rts_run(obs=obs, out_dir=out_dir, tag=tag)
            run_jobs.append(run_job)

        # If not submitting to queue, create the above scripts and bypass the following code
        if not args.no_run:

            # Launch all setup scripts
            for i, setup_job in enumerate(setup_jobs):

                # Setup command and submit to queue
                os.chdir(Path(setup_job).parents[0])
                cmd = f"sbatch {setup_job}"
                setup_job_message = check_output(cmd, shell=True)

                # Use RE to extract job id from output string of setup job
                setup_job_ID = re.search(
                    r"\d+", setup_job_message.decode("utf-8")
                ).group(0)

                # Launch run script with a dependency on the setup scripts
                os.chdir(Path(run_jobs[i]).parents[0])
                cmd = f"sbatch --dependency=afterok:{setup_job_ID} {run_jobs[i]}"
                run_job_message = check_output(cmd, shell=True)
