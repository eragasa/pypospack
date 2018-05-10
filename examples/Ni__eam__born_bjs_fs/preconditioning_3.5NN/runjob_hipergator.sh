#!/bin/sh
#SBATCH --job-name=Ni_fs              # Job name
#SBATCH --mail-type=NONE              # Mail events (NONE, BEGIN, END, FAIL)
#SBATCH --ntasks=128                   # Number of MPI ranks
#SBATCH --cpus-per-task=1             # Number of cores per MPI rank 
#SBATCH --distribution=cyclic:cyclic  # Distribute tasks cyclically on nodes 
#SBATCH --mem-per-cpu=3000mb          # Memory per processor
#SBATCH --time=24:00:00                # Time limit hrs:min:sec
#SBATCH --output=job.out                   # Standard output and error log
#SBATCH --error=job.err
#SBATCH --qos=phillpot-b                # Queue you are submitting to 
pwd; hostname; date
 
module load intel/2016.0.109 
module load openmpi/1.10.2 

echo PYTHONPATH=$PYTHONPATH
echo python=$(which python)
echo PATH=$PATH

echo "start_time:$(date)"
srun --mpi=pmi2 python mc_iterative_sampler.py
echo "end_time:$(date)"
