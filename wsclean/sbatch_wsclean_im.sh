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
timestamp=2021-01-28_2000
data_dir=/astro/mwaeor/MWA/data
out_dir=ws_im
prefix=uvdump_


cd "$data_dir"/"$obsid"/"$timestamp"
mkdir -p $out_dir


time wsclean -name ./"$out_dir"/uvdump_ -scale 0.75amin -size 2048 2048 \
  -niter 10000 -auto-threshold 0.5 -auto-mask 3 \
  -pol I -multiscale -weight briggs 0  -j 20 -mgain 0.85 \
  -no-update-model-required -abs-mem 100 \
  -grid-with-beam -use-idg -idg-mode cpu -pb-undersampling 4 \
  -mwa-path /astro/mwaeor/achokshi/software/local_python/mwa_pb/data \
  "$data_dir"/"$obsid"/"$timestamp"/ms/uvdump_*.ms
