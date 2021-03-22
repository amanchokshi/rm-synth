#!/bin/bash

for i in 1120300232 1120300352 1120082744 1120082864;
do
    for j in ana_wide fee_wide;
    do
        echo Submitting cuffs.sh "$i" "$j" to slurm queue
        sed "s/obsid/"$i"/g" cuffs_parset.in >> "$j"_"$i".in
        sed -i "s/tag/"$j"/g" "$j"_"$i".in
        sbatch cuffs.sh "$j"_"$i".in
    done;
done;