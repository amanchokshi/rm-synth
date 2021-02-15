#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --time=00:20:00
#SBATCH --account=mwaeor
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/uv2ms-%A.out
#SBATCH --error=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/uv2ms-%A.err


obsid=1061316296
tag=patch
data_dir=/astro/mwaeor/achokshi/rm-synth/data
prefix=uvdump_

uv_dir="$data_dir"/"$obsid"/"$tag"/uvfits
ms_dir="$data_dir"/"$obsid"/"$tag"/ms


module use /pawsey/mwa/software/python3/modulefiles
module load casa/5.6.1-8
module load cotter/v4.3
module load intel-mkl/19.0.5

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/pawsey/intel/19.0.5/mkl/lib/intel64/

cd "$uv_dir"

time casa --nologger -c /astro/mwaeor/achokshi/rm-synth/scripts/uv2ms.py $prefix "$obsid"_metafits_ppds.fits

time for band in $(seq -f "%02g" 1 24);
do
  fixmwams "$prefix"${band}.ms "$obsid"_metafits_ppds.fits
done

mkdir -p $ms_dir && mv "$uv_dir"/*.ms $ms_dir
