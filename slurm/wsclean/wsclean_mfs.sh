#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --time=12:00:00
#SBATCH --account=mwaeor
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/wsclean-rm-%A.out
#SBATCH --error=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/wsclean-rm-%A.err

module load wsclean

obsid=1120300352
tag=ana_leakage
data_dir=/astro/mwaeor/achokshi/rm-synth/data
out_dir=wsc_mfs
prefix=uvdump_


cd "$data_dir"/"$obsid"/"$tag"
mkdir -p $out_dir


# MFS
time wsclean -pol iquv -join-polarizations -join-channels \
  -squared-channel-joining --channels-out 24 -weight briggs -1 \
  -scale 0.75amin -size 2048 2048 -mgain 0.95 --niter 10000 \
  -auto-threshold 1 -auto-mask 5 \
  -name ./"$out_dir"/uvdump \
  ./ms/"$prefix"*.ms



# time wsclean -pol iquv -join-polarizations -join-channels \
  # -squared-channel-joining --channels-out 24 -weight briggs -1 \
  # -scale 0.75amin -size 1024 1024 \
  # -multiscale -mgain 0.95 --niter 10000 \
  # -name ./"$out_dir"/uvdump \
  # ./ms/"$prefix"*.ms
  # -no-update-model-required -abs-mem 128 \
  # -grid-with-beam -use-idg -idg-mode cpu -pb-undersampling 4 \
  # -mwa-path /astro/mwaeor/achokshi/software/local_python/mwa_pb/data \
  # -auto-mask 5 -auto-threshold 1 -scale 0.75amin -size 1024 1024 \

# Stokes
# time wsclean -pol QU -join-polarizations -join-channels \
  # -squared-channel-joining --channels-out 768 \
  # -weight briggs -1 \
  # -name ./"$out_dir"/uvdump -scale 0.75amin -size 256 256 \
  # ./ms/"$prefix"*.ms

# Image
# time wsclean -name ./"$out_dir"/uvdump -scale 0.75amin -size 256 256 \
  # -niter 10000 -auto-threshold 0.5 -auto-mask 3 \
  # -pol I -multiscale -weight briggs 0  -j 20 -mgain 0.85 \
  # -no-update-model-required -abs-mem 100 \
  # -grid-with-beam -use-idg -idg-mode cpu -pb-undersampling 4 \
  # -mwa-path /astro/mwaeor/achokshi/software/local_python/mwa_pb/data \
  # "$data_dir"/"$obsid"/"$tag"/ms/"$prefix"*.ms
