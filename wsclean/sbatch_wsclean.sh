#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --time=03:00:00
#SBATCH --account=mwaeor
#SBATCH --nodes=1
#SBATCH --mem=100gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/astro/mwaeor/MWA/data/1061316296/2021-01-20_2000/wsclean-%A.out
#SBATCH --error=/astro/mwaeor/MWA/data/1061316296/2021-01-20_2000/wsclean-%A.err

module load wsclean

cd /astro/mwaeor/MWA/data/1061316296/2021-01-20_2000

mkdir -p images_ii

# time wsclean -name ./images_i/uvdump -size 1024 1024 -niter 10000 \
  # -auto-threshold 0.5 -auto-mask 3 \
  # -pol I -multiscale -weight briggs 0 -scale 0.033 -j 20 -mgain 0.85 \
  # -no-update-model-required \
  # -abs-mem 100 \
  # -grid-with-beam -use-idg -idg-mode cpu -pb-undersampling 4 \
  # -mwa-path /astro/mwaeor/achokshi/software/local_python/mwa_pb/data \
  # /astro/mwaeor/MWA/data/1061316296/2021-01-20_2000/ms/uvdump_*.ms
  # #-channels-out 4 -join-channels \

time wsclean -name ./images_ii/uvdump -size 1024 1024 -niter 10000 \
  -auto-threshold 0.5 -auto-mask 3 -channels-out 24 \
  -pol I -multiscale -weight briggs 0 -scale 0.033 -j 20 -mgain 0.85 \
  -no-update-model-required \
  -abs-mem 100 \
  -grid-with-beam -use-idg -idg-mode cpu -pb-undersampling 4 \
  -mwa-path /astro/mwaeor/achokshi/software/local_python/mwa_pb/data \
  /astro/mwaeor/MWA/data/1061316296/2021-01-20_2000/ms/uvdump_*.ms
