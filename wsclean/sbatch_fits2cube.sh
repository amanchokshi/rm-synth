#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --time=00:30:00
#SBATCH --account=mwaeor
#SBATCH --nodes=1
#SBATCH --mem=30gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/fits2cube-%A.out
#SBATCH --error=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/fits2cube-%A.err


obsid=1061316296
tag=POGSII-EG-024-FEE
data_dir=/astro/mwaeor/achokshi/rm-synth/data
fits_dir="$data_dir"/"$obsid"/"$tag"/wsclean_rm
prefix=uvdump
suffix=dirty
chans=768
dim=256


module load python
module load numpy
module load astropy

mkdir -p "$data_dir"/"$obsid"/"$tag"/rm_cubes

time python /astro/mwaeor/achokshi/rm-synth/wsclean/python-scripts/fits2cube.py \
    --fits_dir=$fits_dir --prefix=$prefix --suffix=$suffix --pol=Q --chans=$chans --dim=$dim --out_dir="$data_dir"/"$obsid"/"$tag"/rm_cubes


time python /astro/mwaeor/achokshi/rm-synth/wsclean/python-scripts/fits2cube.py \
    --fits_dir=$fits_dir --prefix=$prefix --suffix=$suffix --pol=U --chans=$chans --dim=$dim --out_dir="$data_dir"/"$obsid"/"$tag"/rm_cubes
