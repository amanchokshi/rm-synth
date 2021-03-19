
#!/bin/bash

for i in 1120082744 1120082864;
do
    for j in ana_leakage fee_leakage;
    do
        echo Submitting wsclean_mfs.sh "$i" "$j" to slurm queue
        sbatch wsclean_mfs.sh $i $j
    done;
done;
