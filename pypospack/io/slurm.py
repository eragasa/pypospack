""" module to manage jobs to slurm """
import copy
from collections import OrderedDict

slurm_default_configuration_list = [
    ('job_name', 'default_job'),
    ('qos', 'phillpot-b'),
    ('email', 'eragasa@ufl.edu'),
    ('ntasks', '16'),
    ('output', 'job.out'),
    ('error', 'job.err'),
    ('time', '1:00:00'),
    ('memory', '4gb')
]

class SlurmSubmissionScript(object):

    def __init__(self,slurm_dict=None):
        if slurm_dict is None:
            self.configuration = OrderedDict()
        else:
            assert isinstance(slurm_dict,dict)
            self.process_configuration_dictionary(slurm_dict)

    def process_configuration_dictionary(self,slurm_dict):
        assert type(slurm_dict) == OrderedDict
        self.configuration = copy.deepcopy(slurm_dict)

        if 'job_name' not in self.configuration:
            self.configuration['job_name'] = 'default_job'
        if 'qos' not in self.configuration:
            raise ValueError('no queue provided')
        if 'email' not in self.configuration:
            raise ValueError('no email provided')
        if 'ntasks' not in self.configuration:
            self.configuration['ntasks'] = '16'
        if 'output' not in self.configuration:
            self.configuration['output'] = 'job.out'
        if 'error' not in self.configuration:
            self.configuration['error'] = 'job.err'
        if 'time' not in self.configuration:
            self.configuration['time'] = '1:00:00'
        if 'memory' not in self.configuration:
            self.configuration['memory'] = '4gb'

    def section_header_section_to_str(self):
        job_name_ = self.configuration['job_name']
        qos_ = self.configuration['qos']
        mail_ = self.configuration['mail']
        mail_type_ = self.configuration['mail_type']
        ntasks_ = self.configuration['ntasks']
        time_ = self.configuration['time']
        output_ = self.configuration['output']
        error_ = self.configuration['error']
        memory_ = self.configuration['memory']

        str_out = '#!/bin/bash\n'
        str_out += '#SBATCH --job-name={}\n'.format(job_name_)
        str_out += '#SBATCH --qos={}\n'.format(qos_)
        str_out += '#SBATCH --mail-type={}\n'.format(mail_type_)
        str_out += '#SBATCH --mail-user={}\n'.format(mail_)
        str_out += '#SBATCH --ntasks={}\n'.format(ntasks_)
        str_out += '#SBATCH --distribution=cyclic:cyclic\n'
        str_out += '#SBATCH --mem={}\n'.format(memory_)
        str_out += '#SBATCH --time={}\n'.format(time_)
        str_out += '#SBATCH --output={}\n'.format(output_)
        str_out += '#SBATCH --error={}\n'.format(error_)

        return str_out

    def section_load_modules_to_str(self, modules=None):

        # since modules argument is specified, replace modules entry
        if modules is not None:
            self.configuration['modules'] = list(modules)
        modules_ = self.configuration['modules']

        str_out = ""
        for module in modules_:
            str_out += "module load {}\n".format(module)

        return str_out

    def section_slurm_debug_to_str(self):
        str_out = "\n".join([
            'echo slurm_job_id:$SLURM_JOB_ID'
            'echo slurm_job_name:$SLURM_JOB_NAME'
            'echo slurm_job_nodelist:$SLURM_JOB_NODELIST'
            'echo slurm_job_num_nodes:$SLURM_JOB_NUM_NODES'
            'echo slurm_cpus_on_node:$SLURM_CPUS_ON_NODE'
            'echo slurm_ntasks:$SLURM_NTASKS'
            'echo working directory:$(pwd)'
            'echo hostname:$(hostname)'
            'echo start_time:$(date)'
        ]) + "\n"
        return str_out

    def section_command_to_str(self,command=None):
        if command is not None:
            self.configuration['command'] = command
        command_ = self.configuration['command']

        if isinstance(command_, str):
            str_out = command_
        elif isinstance(command_, list):
            str_out = "\n".join(command_) + "\n"
        else:
            raise TypeError()

        return str_out

    def section_postprocessing_to_str(self):
        str_out = "touch jobComplete\n"
        str_out += "echo end_time:$(date)\n"
        return str_out

    def slurm_script_string(self):
        str_out = self.section_header_section_to_str()
        str_out += "#<---------- SLURM debug information"
        str_out += self.section_slurm_debug_to_str()
        str_out += '#<---------- load necessary modules\n'
        str_out += '#<---------- run application\n'
        str_out += self.section_command_to_str()
        str_out += '#<---------- post-processing steps'
        str_out += section_postprocessing_to_str()
        return str_out

    def write(self,filename='runjob.slurm',job_name=None):
        self.filename=filename
        if job_name is not None:
            self.configuration['job_name'] = job_name

        str_out = self.slurm_script_string(self)
        with open(filename,'w') as f:
            f.write(str_out)

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
        error='job.err',
        vasp_bin=None):

    intel_compiler_string = ""
    mpi_compiler_string = ""

    if vasp_bin is None:
        _vasp_bin = os.environ['VASP_BIN']
    else:
        _vasp_bin = vasp_bin

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

    if vasp_bin is not None:
        s += "VASP_BIN={}\n".format(_vasp_bin)
        s += "export VASP_BIN"

    s += 'echo slurm_job_id:$SLURM_JOB_ID\n'
    s += 'echo slurm_job_name:$SLURM_JOB_NAME\n'
    s += 'echo slurm_job_nodelist:$SLURM_JOB_NODELIST\n'
    s += 'echo slurm_job_num_nodes:$SLURM_JOB_NUM_NODES\n'
    s += 'echo slurm_cpus_on_node:$SLURM_CPUS_ON_NODE\n'
    s += 'echo slurm_ntasks:$SLURM_NTASKS\n'

    s += 'echo working directory:$(pwd)\n'
    s += 'echo hostname:$(hostname)\n'
    s += 'echo start_time:$(date)\n'

    s += 'module load intel/2016.0.109\n'
    s += 'module load impi\n'
    s += 'srun --mpi=pmi2 $VASP_BIN > vasp.log\n'
    s += 'echo end_time:$(date)\n'

    with open(filename,'w') as f:
        f.write(s)
