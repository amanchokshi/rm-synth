#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --time=06:00:00
#SBATCH --account=mwaeor
#SBATCH --nodes=1
#SBATCH --mem=128gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/wsclean-rm-%A.out
#SBATCH --error=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/wsclean-rm-%A.err

module load wsclean

obsid=1120300352
tag=fee_leakage
data_dir=/astro/mwaeor/achokshi/rm-synth/data
out_dir=wsc_mfs
prefix=uvdump_


cd "$data_dir"/"$obsid"/"$tag"
mkdir -p $out_dir


time wsclean -pol iquv -join-polarizations -join-channels \
  -squared-channel-joining --channels-out 768 -weight briggs -1 \
  -auto-mask 5 -auto-threshold 1 -scale 0.75amin -size 2048 2048 \
  -multiscale -mgain 0.85 --niter 100000 -no-update-model-required -abs-mem 128 \
  -grid-with-beam -use-idg -idg-mode cpu -pb-undersampling 4 \
  -mwa-path /astro/mwaeor/achokshi/software/local_python/mwa_pb/data \
  -name ./"$out_dir"/uvdump \
  ./ms/"$prefix"*.ms
