#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --time=03:00:00
#SBATCH --account=mwaeor
#SBATCH --nodes=1
#SBATCH --mem=100gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/astro/mwaeor/MWA/data/1061316296/2021-01-27_1200/wsclean-%A.out
#SBATCH --error=/astro/mwaeor/MWA/data/1061316296/2021-01-27_1200/wsclean-%A.err

module load wsclean

obsid=1061316296
timestamp=2021-01-27_1200
data_dir=/astro/mwaeor/MWA/data
out_dir=images_clean
prefix=uvdump_


cd "$data_dir"/"$obsid"/"$timestamp"
mkdir -p $out_dir


time wsclean -pol QU -join-polarizations -join-channels \
  -squared-channel-joining --channels-out 768 \
  -name ./"$out_dir"/uvdump -scale 0.033 -size 1024 1024 \
  "$data_dir"/"$obsid"/"$timestamp"/ms/uvdump_*.ms
