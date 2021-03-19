#!/bin/bash --login

#SBATCH --job-name=wsc_st_2048
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

obsid=1120300352
tag=fee_leakage
data_dir=/astro/mwaeor/achokshi/rm-synth/data
out_dir=wsc_stokes_2048
prefix=uvdump_


cd "$data_dir"/"$obsid"/"$tag"
mkdir -p $out_dir

time wsclean -pol IQUV -join-polarizations -join-channels \
  -squared-channel-joining --channels-out 24 \
  -weight briggs -1 \
  -name ./"$out_dir"/uvdump -scale 0.75amin -size 2048 2048 \
  ./ms/"$prefix"*.ms

# time wsclean -pol QU -join-polarizations -join-channels \
  # -squared-channel-joining --channels-out 768 \
  # -weight briggs -1 \
  # -name ./"$out_dir"/uvdump -scale 0.75amin -size 256 256 \
  # ./ms/"$prefix"*.ms
