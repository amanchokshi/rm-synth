#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --time=04:00:00
#SBATCH --account=mwaeor
#SBATCH --nodes=1
#SBATCH --mem=30gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/noise-%A.out
#SBATCH --error=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/noise-%A.err


obsid=$1
tag=$2
phi_mask=20
cuffs_prefix="$tag"_"$obsid"_
cube_dir=/astro/mwaeor/achokshi/rm-synth/data/"$obsid"/"$tag"/imgs/cubes

module purge
module use /pawsey/mwa/software/python3/modulefiles
module load python-singularity

time python /astro/mwaeor/achokshi/rm-synth/scripts/noise_map.py \
    --phi_mask="$phi_mask" --cuffs_prefix="$cuffs_prefix" --cube_dir="$cube_dir"
