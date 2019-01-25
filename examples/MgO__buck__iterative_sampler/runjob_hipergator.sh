#!/bin/sh
#SBATCH --job-name=Ni_eam             # Job name
#SBATCH --mail-type=NONE              # Mail events (NONE, BEGIN, END, FAIL)
#SBATCH --ntasks=8                    # Number of MPI ranks
#SBATCH --cpus-per-task=1             # Number of cores per MPI rank 
#SBATCH --distribution=cyclic:cyclic  # Distribute tasks cyclically on nodes 
#SBATCH --mem-per-cpu=3000mb          # Memory per processor
#SBATCH --time=8:00:00                # Time limit hrs:min:sec
#SBATCH --output=job.out                   # Standard output and error log
#SBATCH --error=job.err
#SBATCH --qos=phillpot                # Queue you are submitting to 
pwd; hostname; date
 
module load intel/2018.1.163
module load openmpi/3.1.2 

echo PYTHONPATH=$PYTHONPATH
echo PYTHON_BIN=$(which python)
echo PATH=$PATH

export OMPI_MCA_pml=^ucx
export OMPI_MCA_mpi_warn_on_fork=0

srun --mpi=pmi2 $PYTHON_BIN mc_iterative_sampler.py > log.out 
