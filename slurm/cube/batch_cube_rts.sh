#!/bin/bash

# for i in 1120300232 1120300352 1120082744 1120082864;
for i in 1120082744;
do
    for j in 1120082744_+10 1120082744_0 1120082744_-10;
    do
        echo Submitting cube_rts.sh "$i" "$j" to slurm queue
        sbatch cube_rts.sh $i $j
    done;
done;
