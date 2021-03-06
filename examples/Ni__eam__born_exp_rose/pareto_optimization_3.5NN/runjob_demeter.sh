#!/bin/bash
#$ -N pyp_Ni
#$ -cwd
#$ -pe mpi 8
#$ -S /bin/bash
#$ -q single.q
#$ -e job.err
#$ -o job.out

echo start_time:$(date)

# OpenMPI setup
export PATH="/usr/lib64/openmpi/bin:$PATH"
export LD_LIBRARY_PATH="/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH"
MPIRUN=/usr/lib64/openmpi/bin/mpirun

echo The location of mpi is $MPIRUN
echo The python is $(which python)
echo The pythonpath is $PYTHONPATH


/usr/lib64/openmpi/bin/mpirun -np $NSLOTS python dev__eam_Ni__morse_exp_fs.py > log.out

echo stop_time:$(date)
