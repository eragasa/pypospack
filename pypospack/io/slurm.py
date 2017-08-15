""" module to manage jobs to slurm """

def write_vasp_batch_script(filename,job_name,email,qos,ntasks,time,
        output='job.out',
        error='job.err'):
    s = '#!/bin/bash\n'
    s += '#SBATCH --job-name={}\n'.format(job_name)
    s += '#SBATCH --qos={}\n'.format(qos)
    s += '#SBATCH --mail-type=END\n'
    s += '#SBATCH --mail-user={}\n'.format(email)
    s += '#SBATCH --ntasks={}\n'.format(ntasks)
    s += '#SBATCH --distribution=cyclic:cyclic\n'
    s += '#SBATCH --time={}\n'.format(time)
    s += '#SBATCH --output={}\n'.format(output)
    s += '#SBATCH --error={}\n'.format(error)


    s += 'echo slurm_job_id:$SLURM_JOB_ID\n'
    s += 'echo slurm_job_name:$SLURM_JOB_NAME\n'
    s += 'echo slurm_job_nodelist:$SLURM_JOB_NODELIST\n'
    s += 'echo slurm_job_num_nodes:$SLURM_JOB_NUM_NODES\n'
    s += 'echo slurm_cpus_on_node:$SLURM_CPUS_ON_NODE\n'
    s += 'echo slurm_ntasks:$SLURM_NTASKS\n'
    
    s += 'echo working directory:$(pwd)\n'
    s += 'echo hostname:$(hostname)\n'
    s += 'echo start_time:$(date)\n'
    
    s += 'module load intel openmpi vasp\n'
    
    s += 'srun --mpi=pmi2 vasp_std > vasp.log\n'
    
    s += 'echo end_time:$(date)\n'

    with open(filename,'w') as f:
        f.write(s)

