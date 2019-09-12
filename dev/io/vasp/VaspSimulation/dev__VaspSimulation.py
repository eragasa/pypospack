import copy
import math
from pypospack.io.vasp.simulation import VaspSimulation

slurm_default_configuration = [
    ('job_name', 'default_job'),
    ('qos', 'phillpot-b'),
    ('email', 'eragasa@ufl.edu'),
    ('ntasks', '16'),
    ('output', 'job.out'),
    ('error', 'job.err'),
    ('time', '1:00:00'),
    ('memory', '4gb')
]

def determine_npar(n_processors):
    assert isintance(n_processors, int)
    npar = math.floor(math.sqrt(n_processors))

simulation = VaspSimulation()
simulation.submission_script = SlurmSubmissionScript(slurm_default_configuration)
print("n_processors:{}".format(simulation.n_processors))
print("npar:{}".format(determine_npar(simulation.n_processors)))
simulation.npar = determine_npar(n_processors)
