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

module load intel/2016.0.109
module load impi/5.1.1 

export VASP_BIN=/home/eragasa/opt/vasp/vasp.5.4.4/bin/vasp_std

srun --mpi=pmi2 $VASP_BIN > vasp.log
