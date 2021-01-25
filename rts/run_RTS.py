#!/usr/bin/env python

import os
import re
from pathlib import Path
from subprocess import check_output
from sys import exit

import numpy as np
from astropy.io import fits


def make_rts_setup(obs=None, outdir=None, time_stamp=None):
    """Create the rts_setup script for the given observation."""

    # Default source list
    srclist_dir = "/pawsey/mwa/software/python3/srclists/master/"
    srclist_file = "srclist_pumav3_EoR0aegean_EoR1pietro+ForA"

    outfile = open(f"{outdir}{obs}/{time_stamp}/rts_setup_{obs}.sh", "w+")

    outfile.write("#!/bin/bash -l\n")
    outfile.write(f'#SBATCH --job-name="se_{obs}"\n')
    outfile.write(f"#SBATCH --output=RTS-setup-{obs}-%A.out\n")
    outfile.write(f"#SBATCH --error=RTS-setup-{obs}-%A.err\n")
    outfile.write("#SBATCH --nodes=1\n")
    outfile.write("#SBATCH --ntasks-per-node=1\n")
    outfile.write("#SBATCH --time=00:20:00\n")
    outfile.write("#SBATCH --clusters=garrawarla\n")
    outfile.write("#SBATCH --partition=gpuq\n")
    outfile.write("#SBATCH --account=mwaeor\n")
    outfile.write("#SBATCH --gres=gpu:1\n")

    outfile.write("module use /pawsey/mwa/software/python3/modulefiles\n")

    outfile.write("module load python\n")
    outfile.write("module load srclists/master\n")
    outfile.write("module load mongoose\n")
    outfile.write("module load numpy/1.18.2\n")
    outfile.write("module load astropy/4.0.1.post1\n")

    outfile.write(
        "export PYTHONPATH=$PYTHONPATH:/astro/mwaeor/achokshi/software/local_python\n"
    )
    outfile.write("module load mwa_pb\n")

    outfile.write("set -eux\n")
    outfile.write(f"cd {outdir}{obs}/{time_stamp}\n")

    outfile.write("# Generate RTS mwaf files.\n")
    outfile.write("reflag-mwaf-files\n")
    outfile.write("# Generate a source list for the patch step.\n")
    outfile.write("srclist_by_beam.py -n 1000 \\\n")
    outfile.write(f"                   --srclist={srclist_dir}{srclist_file}.txt \\\n")
    outfile.write(
        f"                   --metafits={outdir}{obs}/{time_stamp}/{obs}.metafits\n"
    )

    outfile.write("# Generate a source list for the peel step.\n")
    outfile.write("srclist_by_beam.py -n 3000 \\\n")
    outfile.write(f"                   --srclist={srclist_dir}{srclist_file}.txt \\\n")
    outfile.write(
        f"                   --metafits={outdir}{obs}/{time_stamp}/{obs}.metafits \\\n"
    )
    outfile.write("                   --no_patch \\\n")
    outfile.write("                   --cutoff=90\n")

    outfile.write("# Generate the RTS .in files for both patching and peeling.\n")
    outfile.write("rts-in-file-generator patch \\\n")
    outfile.write(f"                      --base-dir {outdir}{obs}/{time_stamp} \\\n")
    outfile.write("                      --fscrunch 2 \\\n")
    outfile.write(
        f"                      --metafits {outdir}{obs}/{time_stamp}/{obs}.metafits \\\n"
    )
    outfile.write(
        f"                      --srclist {srclist_file}_{obs}_patch1000.txt \\\n"
    )
    outfile.write("                      --num-primary-cals 1 \\\n")
    outfile.write("                      > achokshi_rts_patch.in\n")
    outfile.write("rts-in-file-generator peel \\\n")
    outfile.write(f"                      --base-dir {outdir}{obs}/{time_stamp} \\\n")
    outfile.write("                      --fscrunch 2 \\\n")
    outfile.write(
        f"                      --metafits {outdir}{obs}/{time_stamp}/{obs}.metafits \\\n"
    )
    outfile.write(
        f"                      --srclist {srclist_file}_{obs}_peel3000.txt \\\n"
    )
    outfile.write("                      --num-primary-cals 5 \\\n")
    outfile.write("                      --num-cals 1000 \\\n")
    outfile.write("                      --num-peel 1000 \\\n")
    outfile.write("                      > achokshi_rts_peel.in\n")

    outfile.write("# Ensure permissions are sensible.\n")
    outfile.write("find . -user $USER -type d -exec chmod g+rwx,o+rx,o-w {} \;\n")
    outfile.write("find . -user $USER -type f -exec chmod g+rw,o+r,o-w {} \;\n")
    outfile.write("echo gator rts_setup.sh finished successfully.\n")

    return f"{outdir}{obs}/{time_stamp}/rts_setup_{obs}.sh"


def make_rts_run(obs=None, outdir=None, time_stamp=None):
    """Create the rts_setup script for the given observation."""

    outfile = open(f"{outdir}{obs}/{time_stamp}/rts_run_{obs}.sh", "w+")

    outfile.write("#!/bin/bash -l\n")
    outfile.write(f'#SBATCH --job-name="pa_{obs}"\n')
    outfile.write(f"#SBATCH --output=RTS-patch-{obs}-%A.out\n")
    outfile.write(f"#SBATCH --error=RTS-patch-{obs}-%A.err\n")
    outfile.write("#SBATCH --nodes=25\n")
    outfile.write("#SBATCH --ntasks-per-node=1\n")
    outfile.write("#SBATCH --time=00:30:00\n")
    outfile.write("#SBATCH --clusters=garrawarla\n")
    outfile.write("#SBATCH --partition=gpuq\n")
    outfile.write("#SBATCH --account=mwaeor\n")
    outfile.write("#SBATCH --gres=gpu:1\n")

    outfile.write("module use /pawsey/mwa/software/python3/modulefiles\n")
    outfile.write("module load RTS/sla_to_pal\n")
    outfile.write("module load python-singularity\n")
    # outfile.write('module load RTS/master\n')

    outfile.write("set -eux\n")
    outfile.write("command -v rts_gpu\n")
    outfile.write("export UCX_MEMTYPE_CACHE=n\n")
    outfile.write("date\n")
    outfile.write("srun -n 25 --export=ALL rts_gpu achokshi_rts_patch.in\n")
    outfile.write("date\n")
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
    outfile.write("srun -n 25 --export=ALL rts_gpu achokshi_rts_peel.in\n")
    outfile.write("date\n")

    outfile.write("# Ensure permissions are sensible!\n")
    outfile.write("find . -user $USER -type d -exec chmod g+rwx,o+rx,o-w {} \;\n")
    outfile.write("find . -user $USER -type f -exec chmod g+rw,o+r,o-w {} \;\n")

    return f"{outdir}{obs}/{time_stamp}/rts_run_{obs}.sh"


def test_false(option, value):
    """Test if a necessary argument has a value Exit with a warning message if not."""

    if value is None:
        exit(
            f"{option} must be set for this script - please check your input arguments. Exiting now."
        )
    else:
        pass
    return value


def fhd2rts(metafits, fhd_flags):
    """Convert FHD tile flags to RTS tile flags."""

    # 'pairs' relates the metafits' 'Antenna' column to 'Input'
    # The RTS takes 'Input' // 2 as the tile to be flagged
    hdul = fits.open(metafits)
    pairs = dict(x for x in zip(hdul[1].data["Antenna"], hdul[1].data["Input"]))
    rts_flags = []
    for f in fhd_flags:
        rts_flags.append(str(pairs[int(f)] // 2))
    return rts_flags


def make_flagged_channels(obs=None, outdir=None, time_stamp=None):
    """Create the flagged channel txt file."""

    outfile = open(f"{outdir}{obs}/{time_stamp}/flagged_channels.txt", "w+")

    outfile.write("0\n")
    outfile.write("1\n")
    outfile.write("16\n")
    outfile.write("30\n")
    outfile.write("31\n")

    return f"{outdir}{obs}/{time_stamp}/flagged_channels.txt"


if __name__ == "__main__":

    # Parse the input arguments
    import argparse

    parser = argparse.ArgumentParser(
        description="Setup and run jobs to run the RTS on Garrawarla"
    )
    parser.add_argument(
        "--obs_list",
        default=None,
        required=True,
        help="Text file list of all obs IDs to be processed",
    )
    parser.add_argument(
        "--time_stamp",
        default="2021-01-01_0000",
        help="Timestamp when processing begins - creates timestamp subdir within obsid dir",
    )
    parser.add_argument(
        "--no_run",
        default=False,
        action="store_true",
        help="Don't submit the jobs to the queue",
    )
    parser.add_argument("--fhd_flags", default=None, help="Text file of FHD tile flags")
    parser.add_argument(
        "--check",
        default=False,
        action="store_true",
        help="Run a check on previous specified run for uvfits outputs, then exit",
    )
    parser.add_argument(
        "--clean",
        default=False,
        action="store_true",
        help="Clean out RTS directories of everything but metadata and uvfits, then exit",
    )
    parser.add_argument(
        "--outdir",
        default="/astro/mwaeor/MWA/data/",
        help="Base output dir. Default:/astro/mwaeor/MWA/data/",
    )
    args = parser.parse_args()

    outdir = args.outdir
    time_stamp = args.time_stamp

    # Test the observation list file and make an array of obs
    obs_list = test_false("--obs_list", args.obs_list)
    try:
        obs_list = open(obs_list, "r").read().split("\n")
        obs_list = np.array([obs for obs in obs_list if obs != "" and obs != " "])
    except Exception:
        exit(f"Cannot read --obs_list={obs_list}, please check file. Exiting now.")

    # If set, check the output of all jobs for expected files+size and exit
    if args.check:
        for obs in obs_list:

            # check that 24 uvfits exist
            cmd = f"ls  {outdir}{obs}/{time_stamp}/uvdump_*uvfits | wc -l"
            num_uvfits = int(check_output(cmd, shell=True))

            if num_uvfits < 24:
                print(f"{obs} does not have 24 uvfits files")
            else:
                # now check that they surpass an expected size
                cmd = f"du -c  {outdir}{obs}/{time_stamp}/uvdump_*uvfits | grep total"
                size_uvfits, _ = check_output(cmd, shell=True).split()

                # 2.1G is appropriate for 80kHz, 8 second uvfits.
                if int(size_uvfits) < 2100000:
                    print(f"{obs} uvfits are less than 2.1G")

        exit("Check complete. 24 uvfits files are present and look good.")

    # If set, clean out the RTS directories and exit
    if args.clean:
        for obs in obs_list:

            # Remove gpubox files
            cmd = f"rm -f {outdir}{obs}/{time_stamp}/{obs}*gpubox*fits"
            output = check_output(cmd, shell=True)

            # Remove mwaf files
            cmd = f"rm -f {outdir}{obs}/{time_stamp}/RTS_{obs}_*.mwaf"
            output = check_output(cmd, shell=True)
            cmd = "rm -f {outdir}{obs}/{time_stamp}/{obs}_*.mwaf"
            output = check_output(cmd, shell=True)

            # Temp removal
            cmd = f"rm -f {outdir}{obs}/{time_stamp}/RTS-setup-*-60*"
            output = check_output(cmd, shell=True)
            cmd = f"rm -f {outdir}{obs}/{time_stamp}/RTS-patch-*-60*"
            output = check_output(cmd, shell=True)

        exit("Cleaned directories")

    # Write all shell scripts for running
    setup_jobs = []
    run_jobs = []
    for obs in obs_list:
        txt_job = make_flagged_channels(obs=obs, outdir=outdir, time_stamp=time_stamp)
        setup_job = make_rts_setup(obs=obs, outdir=outdir, time_stamp=time_stamp)
        setup_jobs.append(setup_job)
        run_job = make_rts_run(obs=obs, outdir=outdir, time_stamp=time_stamp)
        run_jobs.append(run_job)

        # Create flagged_tile.txt for each observation by grabbing the latest file (if it exists)
        # and joining it with the list from FHD.
        os.chdir(f"{outdir}{obs}")
        matching = [s for s in os.listdir() if "2019-" in s]
        latest_flag = f"{outdir}{obs}/{max(matching)}/flagged_tiles.txt"

        while not os.path.exists(latest_flag):
            if len(matching) > 1:
                matching.remove(max(matching))
                latest_flag = outdir + obs + "/" + max(matching) + "/flagged_tiles.txt"
            else:
                latest_flag = outdir + obs + "/flagged_tiles.txt"
                break

        # Copy latest flag file
        os.popen("cp " + latest_flag + " " + outdir + obs + "/" + time_stamp)

    if args.fhd_flags:
        fhd_file = args.fhd_flags
        with open(fhd_file) as f:
            for line in f:
                line_arr = line.split()
                obs = line_arr[0]
                # Only write file if flags exist
                if len(line_arr) > 1:
                    fhd_flags = line_arr[1:]
                    metafits_file = (
                        outdir + obs + "/" + time_stamp + "/" + obs + ".metafits"
                    )
                    fhd2rts_flags = fhd2rts(metafits_file, fhd_flags)

                    # Find out if the flagged_tile.txt already exists and get information
                    latest_flag = outdir + obs + "/" + time_stamp + "/flagged_tiles.txt"
                    if os.path.exists(latest_flag):
                        with open(latest_flag) as latest_flag_file:
                            rts_flags = latest_flag_file.read().splitlines()
                            full_flags = list(set(fhd2rts_flags + rts_flags))
                    else:
                        full_flags = fhd2rts_flags

                    with open(latest_flag, "w") as final_flag_file:
                        for flag in full_flags:
                            final_flag_file.write("%s\n" % flag)

    # If not submitting to queue, create the above scripts and bypass the following code
    if args.no_run:
        pass
    else:
        # Launch all setup scripts
        for i, setup_job in enumerate(setup_jobs):

            # Setup command and submit to queue
            os.chdir(Path(setup_job).parents[0])
            cmd = f"sbatch  {setup_job}"
            setup_job_message = check_output(cmd, shell=True)

            # Use RE to extract job id from output string of setup job
            setup_job_ID = re.search(r"\d+", setup_job_message.decode('utf-8')).group(0)

            # Launch run script with a dependency on the setup scripts
            os.chdir(Path(run_jobs[i]).parents[0])
            cmd = f"sbatch --dependency=afterok:{setup_job_ID} {run_jobs[i]}"
            run_job_message = check_output(cmd, shell=True)
