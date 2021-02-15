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


fits_dir=/astro/mwaeor/achokshi/rm-synth/iono/tmp_fits
out_dir=/astro/mwaeor/achokshi/rm-synth/iono/tmp_rot
pogs_path=/astro/mwaeor/achokshi/rm-synth/iono/POGS-II_ExGal.fits

# Clean slate
module purge

# The epic singularity image
module use /astro/mwaeor/achokshi/software/modulefiles/
module load rmextract-singularity


time python /astro/mwaeor/achokshi/rm-synth/iono/corr_iono.py \
    --fits_dir=$fits_dir --out_dir=$out_dir --pogs_path=$pogs_path
