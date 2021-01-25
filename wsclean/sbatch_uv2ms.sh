#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --partition=workq
#SBATCH --time=00:20:00
#SBATCH --account=mwaeor
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module use /pawsey/mwa/software/python3/modulefiles
module load casa/5.6.1-8
module load cotter/v4.3
module load intel-mkl/19.0.5

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/pawsey/intel/19.0.5/mkl/lib/intel64/

cd /astro/mwaeor/MWA/data/1061316296/2021-01-20_2000/

time casa --nologger -c /astro/mwaeor/achokshi/rm-synth/wsclean/uv2ms.py uvdump_ /astro/mwaeor/MWA/data/1061316296/2021-01-20_2000/1061316296_metafits_ppds.fits


time for band in '01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24'
do
  fixmwams uvdump_${band}.ms /astro/mwaeor/MWA/data/1061316296/2021-01-20_2000/1061316296_metafits_ppds.fits
done

mkdir ms && mv *.ms ms
