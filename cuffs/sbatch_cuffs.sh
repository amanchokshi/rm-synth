#!/bin/bash -l

#SBATCH --job-name="cuFFS"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --clusters=garrawarla
#SBATCH --partition=gpuq
#SBATCH --account=mwaeor
#SBATCH --gres=gpu:1
#SBATCH --output=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/cuffs-%A.out
#SBATCH --error=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/cuffs-%A.err

module use /pawsey/mwa/software/python3/modulefiles
module load RTS/sla_to_pal

module use /astro/mwaeor/achokshi/software/modulefiles
module load cuFFS

time rmsynthesis /astro/mwaeor/achokshi/rm-synth/cuffs/parsetFile
