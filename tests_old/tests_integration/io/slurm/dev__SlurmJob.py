from collections import OrderedDict
import pypospack.io.slurm as slurm

if __name__ == "__main__":
    
    slurm_dict = OrderedDict()
    slurm_dict['email'] = 'eragasa@ufl.edu'
    slurm_dict['qos'] = 'phillpot'
    slurm_dict['ntasks'] = 16
    slurm_dict['time'] = "1:00:00"
    slurm_dict['output'] = 'job.out'
    slurm_dict['error'] = 'job.err'
    slurm_dict['modules'] = [
        'intel/2018.1.163',
        'openmpi/2.0.3',
        'vasp/5.4.4']
    slurm_dict['command'] = 'vasp_std > vasp.log'
 

    _slurm_filename = 'runjob.slurm'
    slurm_job = slurm.SlurmSubmissionScript(slurm_dict)
    slurm_job.write(filename=_slurm_filename)

