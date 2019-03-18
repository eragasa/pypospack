#!/bin/sh
#SBATCH --job-name=Al_fs              # Job name
#SBATCH --mail-type=NONE              # Mail events (NONE, BEGIN, END, FAIL)
#SBATCH --ntasks=64                   # Number of MPI ranks
#SBATCH --cpus-per-task=1             # Number of cores per MPI rank 
#SBATCH --distribution=cyclic:cyclic  # Distribute tasks cyclically on nodes 
#SBATCH --mem-per-cpu=3000mb          # Memory per processor
#SBATCH --time=12:00:00                # Time limit hrs:min:sec
#SBATCH --output=job.out                   # Standard output and error log
#SBATCH --error=job.err
#SBATCH --qos=phillpot                # Queue you are submitting to 
pwd; hostname; date

module load intel/2016.0.109
module load impi/5.1.1

echo PYTHONPATH=$PYTHONPATH
PYTHON_BIN=$(which python)
echo PYTHON_BIN=$PYTHON_BIN
echo python=$(which python)
echo PATH=$PATH

echo "start_time:$(date)"
#export OMPI_MCA_pml=^ucx
srun --mpi=pmi2 $PYTHON_BIN mc_iterative_sampler.py
# srun python mc_iterative_sampler.py --mpi=pmix_v1
# mpiexec python mc_iterative_sampler.py
echo "end_time:$(date)"

