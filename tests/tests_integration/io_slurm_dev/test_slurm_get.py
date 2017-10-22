import subprocess, shlex, copy

def get_jobs_from_slurm_running_only(username):
    cmd = 'squeue -u {} -t RUNNING'.format(username)
    args = shlex.split(cmd)
    output = subprocess.run(args,stdout=subprocess.PIPE)
    return output.stdout.decode('utf-8'

def get_jobs_from_slurm_pending_only(username):
    cmd = 'squeue -u {} -t PENDING'.format(username)
    args = shlex.split(cmd)
    output = subprocess.run(args,stdout=subprocess.PIPE)
    return output.stdout.decode('utf-8')

def get_jobs_from_slurm(username):
    cmd = 'squeue -u {}'.format(username)
    args = shlex.split(cmd)
    output = subprocess.run(args,stdout=subprocess.PIPE)
    return output.stdout.decode('utf-8')

def get_job_control_info(jobid):
    cmd = 'control show jobid -dd {}'.format(jobid)
    args = shlex.split(cmd)
    output = subprocess.run(args,stdout=subprocess.PIPE)
    return output.stdout.decode('utf-8')

def process_get_jobs_info(results):
    lines = results.split('\n')
    
    names = None # initialize
    jobinfo = [] # initialize
    for i,line in enumerate(lines):
        if i == 0:
            names = [s.strip() for s in line.split()]
        else:
            result.append([s.strip() for s in line.split()])
    return_dict = {
        'names':names,
        'jobinfo':jobinfo
        }
    return copy.deepcopy(return_dict)

username = 'eragasa'
output = get_jobs_from_slurm(username)
