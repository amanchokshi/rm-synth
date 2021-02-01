#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --time=01:00:00
#SBATCH --account=mwaeor
#SBATCH --nodes=1
#SBATCH --mem=100gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/wsclean-im-%A.out
#SBATCH --error=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/wsclean-im-%A.err

module load wsclean

obsid=1061316296
tag=POGSII-EG-024
data_dir=/astro/mwaeor/achokshi/rm-synth/data
out_dir=wsclean_im
prefix=uvdump_


cd "$data_dir"/"$obsid"/"$tag"
mkdir -p $out_dir


time wsclean -name ./"$out_dir"/uvdump -scale 0.75amin -size 256 256 \
  -niter 10000 -auto-threshold 0.5 -auto-mask 3 \
  -pol I -multiscale -weight briggs 0  -j 20 -mgain 0.85 \
  -no-update-model-required -abs-mem 100 \
  -grid-with-beam -use-idg -idg-mode cpu -pb-undersampling 4 \
  -mwa-path /astro/mwaeor/achokshi/software/local_python/mwa_pb/data \
  "$data_dir"/"$obsid"/"$tag"/ms/"$prefix"*.ms
