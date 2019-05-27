#!/bin/bash
#
#$ -N MgO
#$ -cwd
#$ -S /bin/bash
#$ -pe mpi 8 
#$ -q all.q,eight.q,four.q,single.q
#@compute-0-11,all.q@compute-0-0,all.q@compute-0-2,all.q@compute-0-25,all.q@compute-0-26,all.q@compute-0-28


#unset SGE_ROOT

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/mkl/10.0.011/lib/em64t:/opt/intel/fc/9.1.040/lib:$HOME/lib/32:$HOME/lib

export MKL_NUM_THREADS=1

#$ -v $LD_LIBRARY_PATH
#$ -v $MKL_NUM_THREADS
export RUNPATH=`pwd`

/opt/intel/mpich2-1.4.1p1/bin/mpiexec -rmk pbs -machinefile $TMPDIR/machines -np $NSLOTS $RUNPATH/PhonTS > out.dat

#/opt/intel/mpich2-1.4.1p1/bin/mpiexec -rmk pbs -machinefile $TMPDIR/machines -np $NSLOTS ~/SOFT/PhonTS-1.0.3/src/PhonTS > out.dat



