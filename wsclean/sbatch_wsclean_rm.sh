#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --time=01:00:00
#SBATCH --account=mwaeor
#SBATCH --nodes=1
#SBATCH --mem=100gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/wsclean-rm-%A.out
#SBATCH --error=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/wsclean-rm-%A.err

module load wsclean

obsid=1061316296
tag=POGSII-EG-024
data_dir=/astro/mwaeor/achokshi/rm-synth/data
out_dir=wsclean_rm
prefix=uvdump_


cd "$data_dir"/"$obsid"/"$tag"
mkdir -p $out_dir


time wsclean -pol QU -join-polarizations -join-channels \
  -squared-channel-joining --channels-out 768 \
  -name ./"$out_dir"/uvdump -scale 0.75amin -size 256 256 \
  ./ms/"$prefix"*.ms
