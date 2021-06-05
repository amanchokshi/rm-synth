#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --time=1-00
#SBATCH --account=mwaeor
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/beam_mcmc-%A.out
#SBATCH --error=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/beam_mcmc-%A.err

module load python
module load numpy
module load healpy

export PYTHONPATH=$PYTHONPATH:/astro/mwaeor/achokshi/software/local_python/
export MWA_BEAM_FILE=/astro/mwaeor/achokshi/software/local_python/mwa_pb/data/mwa_full_embedded_element_pattern.h5

module load hyperbeam
module load emcee

time python /astro/mwaeor/achokshi/rm-synth/scripts/beam_mcmc.py \
    --map_dir="/astro/mwaeor/achokshi/rm-synth/data/embers_maps" \
    --map_name="S06YY_rf1YY" --out_dir="/astro/mwaeor/achokshi/rm-synth/data/beam_mcmc"
