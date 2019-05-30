#!/bin/bash

#SBATCH -J Argon_test
#SBATCH --ntasks=16
#SBATCH --out=Forge-%j.out
#SBATCH --time=0-00:20:00
#SBATCH --mail-type=begin,end,fail,requeue

module load intel/2015.2.164
module load mvapich2/intel/15/eth
export RUNPATH=`pwd`
mpirun $RUNPATH/PhonTS > out.dat


