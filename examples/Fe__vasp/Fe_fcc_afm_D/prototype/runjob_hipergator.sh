#!/bin/sh
#SBATCH --job-name=Fe_fcc_vacA         # Job name
#SBATCH --mail-type=NONE              # Mail events (NONE, BEGIN, END, FAIL)
#SBATCH --ntasks=64                   # Number of MPI ranks
#SBATCH --cpus-per-task=1             # Number of cores per MPI rank 
#SBATCH --distribution=cyclic:cyclic  # Distribute tasks cyclically on nodes 
#SBATCH --time=8:00:00                # Time limit hrs:min:sec
#SBATCH --output=job.out              # Standard output and error log
#SBATCH --error=job.err
#SBATCH --qos=phillpot-b              # Queue you are submitting to 
pwd; hostname; date

module load intel/2019.1.144  
module load openmpi/4.0.0
module load vasp/5.4.4

srun --mpi=pmi2 vasp_std > vasp.log
