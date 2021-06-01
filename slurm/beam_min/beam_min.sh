#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --time=12:00:00
#SBATCH --account=mwaeor
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/beam_min-%A.out
#SBATCH --error=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/beam_min-%A.err

module load python
module load numpy
module load scipy
module load healpy

export PYTHONPATH=$PYTHONPATH:/astro/mwaeor/achokshi/software/local_python/
export MWA_BEAM_FILE=/astro/mwaeor/achokshi/software/local_python/mwa_pb/data/mwa_full_embedded_element_pattern.h5

module load hyperbeam

time python /astro/mwaeor/achokshi/rm-synth/scripts/beam_minimize.py \
    --sat_map="$1"
