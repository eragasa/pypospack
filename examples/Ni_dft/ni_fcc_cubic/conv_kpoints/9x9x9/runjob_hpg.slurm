#!/bin/bash
#SBATCH --job-name=Ni_fcc_9x9x9
#SBATCH --qos=phillpot
#SBATCH --mail-type=END
#SBATCH --mail-user=sullberg@ufl.edu
#SBATCH --ntasks=16
#SBATCH --distribution=cyclic:cyclic
#SBATCH --time=1:00:00
#SBATCH --output=job.out
#SBATCH --error=job.err
echo slurm_job_id:$SLURM_JOB_ID
echo slurm_job_name:$SLURM_JOB_NAME
echo slurm_job_nodelist:$SLURM_JOB_NODELIST
echo slurm_job_num_nodes:$SLURM_JOB_NUM_NODES
echo slurm_cpus_on_node:$SLURM_CPUS_ON_NODE
echo slurm_ntasks:$SLURM_NTASKS
echo working directory:$(pwd)
echo hostname:$(hostname)
echo start_time:$(date)
module load intel/2016.0.109
module load impi
srun --mpi=pmi2 $VASP_BIN > vasp.log
echo end_time:$(date)
