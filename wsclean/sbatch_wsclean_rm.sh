#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --time=03:00:00
#SBATCH --account=mwaeor
#SBATCH --nodes=1
#SBATCH --mem=100gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/astro/mwaeor/MWA/data/1061316296/2021-01-28_2000/wsclean-%A.out
#SBATCH --error=/astro/mwaeor/MWA/data/1061316296/2021-01-28_2000/wsclean-%A.err

module load wsclean

obsid=1061316296
tag=2021-01-31_0100
data_dir=/astro/mwaeor/achokshi/rm-synth/data
out_dir=ws_rm
prefix=uvdump_


cd "$data_dir"/"$obsid"/"$tag"
mkdir -p $out_dir


time wsclean -pol QU -join-polarizations -join-channels \
  -squared-channel-joining --channels-out 768 \
  -name ./"$out_dir"/uvdump -scale 0.75amin -size 2048 2048 \
  "$data_dir"/"$obsid"/"$timestamp"/ms/uvdump_*.ms
