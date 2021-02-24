#!/bin/bash -l
#SBATCH --job-name="download"
#SBATCH --export=None
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --mem=1G
#SBATCH --clusters=garrawarla
#SBATCH --partition=workq
#SBATCH --account=mwaeor
#SBATCH --export=MWA_ASVO_API_KEY
#SBATCH --output=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/download-%A.out
#SBATCH --error=/astro/mwaeor/achokshi/rm-synth/data/slurm-logs/download-%A.err

module use /pawsey/mwa/software/python3/modulefiles
module load giant-squid

module use /astro/mwaeor/achokshi/software/modulefiles
module load jq

sub_dir=gpubox
mwa_dir=/astro/mwaeor/MWA/data
obs_file=/astro/mwaeor/achokshi/rm-synth/slurm/rts/rts_in.json

# Read obsids from json
obs_id_array=$(jq '.obsids | .[]' "$obs_file")

for obs_id in "${obs_id_array[@]}"
do
   giant-squid submit-vis "$obs_id"
done


for obs_id in "${obs_id_array[@]}"
do
    mkdir -p "${mwa_dir}"/"${obs_id}"/"${sub_dir}"
    cd "${mwa_dir}"/"${obs_id}"/"${sub_dir}"
    echo Made and moved to "${mwa_dir}"/"${obs_id}"/"${sub_dir}"

    num_files=$(ls "${mwa_dir}"/"${obs_id}"/"${sub_dir}" -1 | wc -l)

    if [[ $num_files -gt 47 ]]; then
        echo "$obs_id" was already downloaded
    else
        download_flag=0
        while [ $download_flag -eq 0 ]
        do
            message=$(giant-squid download "${obs_id}")

            if [[ $message == *"not ready"* ]]; then
                echo "$(date) - "$num_files" downloaded"
                sleep 60
            else
                unzip "${mwa_dir}"/"${obs_id}"/"${sub_dir}"/"${obs_id}"_flags.zip
                rm "${mwa_dir}"/"${obs_id}"/"${sub_dir}"/"${obs_id}"_flags.zip
                download_flag=1
                echo Successfully downloaded "${obs_id}"
            fi
        done
    fi
done
