#!/bin/bash

for i in 1120300232 1120300352 1120082744 1120082864;
do
    for j in ana_wide fee_wide;
    do
        echo Submitting cube_rts.sh "$i" "$j" to slurm queue
        sbatch cube_rts.sh $i $j
    done;
done;
