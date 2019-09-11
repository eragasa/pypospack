import pytest
from collections import OrderedDict
from pypospack.io.slurm import SlurmSubmissionScript

slurm_configuration_list = [
    ['job_name', 'job_name'],
    ['qos', 'phillpot'],
    ['email', 'eragasa@ufl.edu'],
    ['ntasks', 32],
    ['output', 'job.out'],
    ['error', 'job.err'],
    ['time','1:00:00'],
    ['memory', '4gb']
]

slurm_configuration_dict = OrderedDict(
    [(k[0],k[1]) for k in slurm_configuration_list]
)

if __name__ == "__main__":
    slurm_script = SlurmSubmissionScript()
    assert isinstance(slurm_script, SlurmSubmissionScript)
    assert isinstance(slurm_script.configuration, OrderedDict)

    slurm_script.configuration = slurm_configuration_dict
    assert isinstance(slurm_script.configuration, OrderedDict)
    slurm_header_section = slurm_script.section_header_section_to_str()
    print(slurm_header_section)
