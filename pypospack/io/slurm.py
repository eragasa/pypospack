""" module to manage jobs to slurm """
import copy
from collections import OrderedDict

class SlurmSubmissionScript(object):

    def __init__(self,slurm_dict=None):
        if slurm_dict is None:
            self.configuration = OrderedDict()
        else:
            assert isinstance(slurm_dict,dict)
            self.process_configuration_dictionary(slurm_dict)

    def process_configuration_dictionary(self,slurm_dict):
        assert type(slurm_dict) == OrderedDict
        self.configuration =copy.deepcopy(slurm_dict)

        if 'job_name' not in self.configuration:
            self.configuration['job_name'] = 'default_job'
        if 'qos' not in self.configuration:
            raise ValueError('no queue provided')
        if 'email' not in self.configuration:
            raise ValueError('no email provided')
        if 'ntasks' not in self.configuration:
            self.configuration = '16'
        if 'output' not in self.configuration:
            self.configuration['output'] = 'job.out'
        if 'error' not in self.configuration:
            self.configuration['error'] = 'job.err'
        if 'time' not in self.configuration:
            self.configuration['time'] = '1:00:00'

    def write(self,filename='runjob.slurm',job_name='default'):
        self.filename=filename
       
        self.configuration['job_name'] = job_name 
        _job_name = self.configuration['job_name']
        _qos = self.configuration['qos']
        _email = self.configuration['email']
        _ntasks = self.configuration['ntasks']
        _time = self.configuration['time']
        _output = self.configuration['output']
        _error = self.configuration['error']
        _cmd = self.configuration['command']

        s = '#!/bin/bash\n'
        s += '#SBATCH --job-name={}\n'.format(_job_name)
        s += '#SBATCH --qos={}\n'.format(_qos)
        s += '#SBATCH --mail-type=END\n'
        s += '#SBATCH --mail-user={}\n'.format(_email)
        s += '#SBATCH --ntasks={}\n'.format(_ntasks)
        s += '#SBATCH --distribution=cyclic:cyclic\n'
        s += '#SBATCH --time={}\n'.format(_time)
        s += '#SBATCH --output={}\n'.format(_output)
        s += '#SBATCH --error={}\n'.format(_error)
        s += 'echo slurm_job_id:$SLURM_JOB_ID\n'
        s += 'echo slurm_job_name:$SLURM_JOB_NAME\n'
        s += 'echo slurm_job_nodelist:$SLURM_JOB_NODELIST\n'
        s += 'echo slurm_job_num_nodes:$SLURM_JOB_NUM_NODES\n'
        s += 'echo slurm_cpus_on_node:$SLURM_CPUS_ON_NODE\n'
        s += 'echo slurm_ntasks:$SLURM_NTASKS\n'
        s += 'echo working directory:$(pwd)\n'
        s += 'echo hostname:$(hostname)\n'
        s += 'echo start_time:$(date)\n'

        s += '#<---------- load necessary modules\n'
        if 'modules' in self.configuration:
            for module_name in self.configuration['modules']:
                s += "module load {}\n".format(module_name)	
	
        s += '#<---------- run application\n'
        s += 'srun --mpi=pmi2 {}\n'.format(_cmd)
        s += 'touch jobCompleted\n'    
        s += 'echo end_time:$(date)\n'

        with open(filename,'w') as f:
            f.write(s)

def write_phonts_batch_script(filename,job_name,email,qos,ntasks,time,
        output='job.out',error='job.err'):
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
    
    s += 'module load intel openmpi\n'
    s += '\n'
    s += 'srun --mpi=pmi2 $PHONTS_BIN > phonts.log\n'
    s += 'touch jobCompleted\n'    
    s += 'echo end_time:$(date)\n'

    with open(filename,'w') as f:
        f.write(s)

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
    s += 'touch jobCompleted\n'    
    s += 'echo end_time:$(date)\n'

    with open(filename,'w') as f:
        f.write(s)

