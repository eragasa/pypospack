#!/bin/sh
#SBATCH --job-name=hello_world_test   # Job name
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
echo python=$(which python)
echo PATH=$PATH

echo "start_time:$(date)"
srun --mpi=pmix_v1 python mpi_code.py
echo "end_time:$(date)"
