#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --time=00:30:00
#SBATCH --account=mwaeor
#SBATCH --nodes=1
#SBATCH --mem=30gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/cube_rts-%A.out
#SBATCH --error=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/cube_rts-%A.err


obsid=1120300352
data_dir=/astro/mwaeor/achokshi/rm-synth/data/"$obsid"/"$1"/imgs

module load python
module load numpy
module load astropy

time python /astro/mwaeor/achokshi/rm-synth/scripts/cube_rts.py \
    --fits_dir="$data_dir"/stokes --out_dir="$data_dir"/cubes

time python /astro/mwaeor/achokshi/rm-synth/scripts/cube_rts.py \
    --fits_dir="$data_dir"/stokes_iono --out_dir="$data_dir"/cubes_iono
