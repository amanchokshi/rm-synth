#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --time=00:20:00
#SBATCH --account=mwaeor
#SBATCH --nodes=1
#SBATCH --mem=30gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/noise-%A.out
#SBATCH --error=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/noise-%A.err


obsid=1120300232
phi_mask=20
cuffs_prefix=rts_imgr_
cube_dir=/astro/mwaeor/achokshi/rm-synth/data/"$obsid"/"$1"/imgs/cubes_iono

module load python
module load numpy
module load scipy
module load astropy
module load matplotlib

time python /astro/mwaeor/achokshi/rm-synth/scripts/noise_map.py \
    --phi_mask="$phi_mask" --cuffs_prefix="$cuffs_prefix" --cube_dir="$cube_dir"
