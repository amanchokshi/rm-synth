#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --time=00:10:00
#SBATCH --account=mwaeor
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/astro/mwaeor/MWA/data/1061316296/2021-01-20_2000/fits2cube-%A.out
#SBATCH --error=/astro/mwaeor/MWA/data/1061316296/2021-01-20_2000/fits2cube-%A.err


obsid=1061316296
timestamp=2021-01-20_2000
data_dir=/astro/mwaeor/MWA/data
fits_dir="$data_dir"/"$obsid"/"$timestamp"/images_iii
prefix=uvdump
chans=384
dim=1024

module load python
module load numpy
module load astropy

export PYTHONPATH=$PYTHONPATH:/astro/mwaeor/achokshi/software/local_python
# module load spectral-cube

time python /astro/mwaeor/achokshi/rm-synth/wsclean/python-scripts/fits2cube.py --fits_dir=$fits_dir --prefix=$prefix --pol=Q --chans=$chans --dim=$dim --out_dir=$fits_dir

time python /astro/mwaeor/achokshi/rm-synth/wsclean/python-scripts/fits2cube.py --fits_dir=$fits_dir --prefix=$prefix --pol=U --chans=$chans --dim=$dim --out_dir=$fits_dir
