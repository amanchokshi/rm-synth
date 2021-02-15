#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --time=00:10:00
#SBATCH --account=mwaeor
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/corr_iono-%A.out
#SBATCH --error=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/corr_iono-%A.err


rts_dir=/astro/mwaeor/achokshi/rm-synth/data/1061316296/rts_imgr/rts_iii
pogs_path=/astro/mwaeor/achokshi/rm-synth/slurm/iono/POGS-II_ExGal.fits

# Clean slate
module purge

# The epic singularity image
module use /astro/mwaeor/achokshi/software/modulefiles/
module load rmextract-singularity


time python /astro/mwaeor/achokshi/rm-synth/scripts/corr_iono.py \
    --fits_dir="$rts_dir"/stokes --out_dir="$rts_dir"/stokes_iono --pogs_path=$pogs_path
