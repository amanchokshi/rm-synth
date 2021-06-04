#!/bin/bash

map_dir="/astro/mwaeor/achokshi/rm-synth/data/embers_maps"
out_dir="/astro/mwaeor/achokshi/rm-synth/data/beam_min_1024_masked"
beam_min="/astro/mwaeor/achokshi/rm-synth/slurm/beam_min/beam_min.sh"

map_arr=( "S06XX_rf1XX" "S06YY_rf1YY" "S07XX_rf1XX" "S07YY_rf1YY" "S08XX_rf1XX" "S08YY_rf1YY" "S09XX_rf1XX" "S09YY_rf1YY" "S10XX_rf1XX" "S10YY_rf1YY" "S12XX_rf1XX" "S12YY_rf1YY" "S29XX_rf1XX" "S29YY_rf1YY" "S30XX_rf1XX" "S30YY_rf1YY" "S31XX_rf1XX" "S31YY_rf1YY" "S32XX_rf1XX" "S32YY_rf1YY" "S33XX_rf1XX" "S33YY_rf1YY" "S34XX_rf1XX" "S34YY_rf1YY" "S35XX_rf1XX" "S35YY_rf1YY" "S36XX_rf1XX" "S36YY_rf1YY" )

for m in ${map_arr[*]};
do
    sbatch ${beam_min} ${map_dir} ${m} ${out_dir}
done
