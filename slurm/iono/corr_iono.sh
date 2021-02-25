#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --time=00:30:00
#SBATCH --account=mwaeor
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/corr_iono-%A.out
#SBATCH --error=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/corr_iono-%A.err

# Input params
obsid=1120300232
pogs_obj=POGSII-EG-321
pogs_path=/astro/mwaeor/achokshi/rm-synth/slurm/iono/POGS-II_ExGal.fits
data_dir=/astro/mwaeor/achokshi/rm-synth/data/"$obsid"/"$1"/imgs

# Clean slate
module purge

# The epic singularity image
module use /astro/mwaeor/achokshi/software/modulefiles/
module load rmextract-singularity


time python /astro/mwaeor/achokshi/rm-synth/scripts/corr_iono.py \
    --obsid=$obsid --fits_dir="$data_dir"/stokes --out_dir="$data_dir"/stokes_iono --pogs_obj=$pogs_obj --pogs_path=$pogs_path
