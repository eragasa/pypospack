#!/bin/bash
#$ -N Si_neb
#$ -cwd
#$ -pe mpi 4
#$ -S /bin/bash
#$ -q all2.q
#$ -q all.q
#$ -q eight.q
#$ -e $JOB_NAME.e$JOB_ID
#$ -o $JOB_NAME.o$JOB_ID

lammps_bin=lmps_bin
mpiexec=mpi
param_file=culled_009.out
n_sim=10
for (( i=0; i<$n_sim; i++))
    do
	python read_data_file.py min $i $param_file
	echo "${mpiexec} $lammps_bin -in in.min"
	python read_data_file.py neb $i $param_file
	echo "${mpiexec} -np 4 $lammps_bin -partition 4x1 -in in.neb"
	python read_data_file.py post $i $param_file
    done
