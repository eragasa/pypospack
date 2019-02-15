#!/bin/bash

for d in */; do
    echo $d
    cp ../minimize_init/CONTCAR $d/POSCAR
    cd $d
    bash vaspcleanup.sh
    sbatch runjob_hpg.slurm
    cd ../
done

