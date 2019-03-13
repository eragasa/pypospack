#!/bin/sh
#SBATCH --job-name=ompi_environment   # Job name
#SBATCH --mail-type=NONE              # Mail events (NONE, BEGIN, END, FAIL)
#SBATCH --ntasks=4                    # Number of MPI ranks
#SBATCH --cpus-per-task=1             # Number of cores per MPI rank 
#SBATCH --distribution=cyclic:cyclic  # Distribute tasks cyclically on nodes 
#SBATCH --mem-per-cpu=3000mb          # Memory per processor
#SBATCH --time=1:00:00                # Time limit hrs:min:sec
#SBATCH --output=job.out              # Standard output and error log
#SBATCH --error=job.err
#SBATCH --qos=phillpot                # Queue you are submitting to 
pwd; hostname; date
 
module load intel/2018.1.163 
module load openmpi/3.1.2

export OMPI_MCA_mpi_warn_on_fork=0

echo PYTHONPATH=$PYTHONPATH

export PYTHON_BIN=$(which python)
echo PYTHON_BIN=$PYTHON_BIN
echo PATH=$PATH

echo "start_time:$(date)"
export OMPI_MCA_pml=^ucx
srun --mpi=pmi2 $PYTHON_BIN test__ompi_environment.py   
echo "end_time:$(date)"
