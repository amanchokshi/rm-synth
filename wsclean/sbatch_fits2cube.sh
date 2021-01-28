#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --time=00:20:00
#SBATCH --account=mwaeor
#SBATCH --nodes=1
#SBATCH --mem=20gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/astro/mwaeor/MWA/data/1061316296/2021-01-28_2000/fits2cube-%A.out
#SBATCH --error=/astro/mwaeor/MWA/data/1061316296/2021-01-28_2000/fits2cube-%A.err


obsid=1061316296
timestamp=2021-01-28_2000
data_dir=/astro/mwaeor/MWA/data
fits_dir="$data_dir"/"$obsid"/"$timestamp"/ws_rm
prefix=uvdump
suffix=dirty
chans=768
dim=2048

module load python
module load numpy
module load astropy

mkdir -p "$data_dir"/"$obsid"/"$timestamp"/rm_cubes

time python /astro/mwaeor/achokshi/rm-synth/wsclean/python-scripts/fits2cube.py \
    --fits_dir=$fits_dir --prefix=$prefix --suffix=$suffix --pol=Q --chans=$chans --dim=$dim --out_dir="$data_dir"/"$obsid"/"$timestamp"/rm_cubes


time python /astro/mwaeor/achokshi/rm-synth/wsclean/python-scripts/fits2cube.py \
    --fits_dir=$fits_dir --prefix=$prefix --suffix=$suffix --pol=U --chans=$chans --dim=$dim --out_dir="$data_dir"/"$obsid"/"$timestamp"/rm_cubes

