import pytest
from collections import OrderedDict
from pypospack.io.slurm import SlurmSubmissionScript

def dev__section_header_to_str__no_args():
    slurm_configuration_list = [
        ['job_name', 'job_name'],
        ['qos', 'phillpot'],
        ['mail_type', 'END'],
        ['mail', 'eragasa@ufl.edu'],
        ['ntasks', 32],
        ['output', 'job.out'],
        ['error', 'job.err'],
        ['time','1:00:00'],
        ['memory', '4gb']
    ]

    slurm_configuration_dict = OrderedDict(
        [(k[0],k[1]) for k in slurm_configuration_list]
    )

    slurm_header_str = "\n".join([
        "#!/bin/bash",
        "#SBATCH --job-name=job_name",
        "#SBATCH --qos=phillpot",
        "#SBATCH --mail-type=END",
        "#SBATCH --mail-user=eragasa@ufl.edu",
        "#SBATCH --ntasks=32",
        "#SBATCH --distribution=cyclic:cyclic",
        "#SBATCH --mem=4gb",
        "#SBATCH --time=1:00:00",
        "#SBATCH --output=job.out",
        "#SBATCH --error=job.err"
    ]) + "\n"

    slurm_script = SlurmSubmissionScript()
    assert isinstance(slurm_script, SlurmSubmissionScript)
    assert isinstance(slurm_script.configuration, OrderedDict)

    slurm_script.configuration = slurm_configuration_dict
    assert isinstance(slurm_script.configuration, OrderedDict)
    slurm_header_section = slurm_script.section_header_section_to_str()
    print(slurm_header_section)
    print(slurm_header_str)
    assert slurm_header_section == slurm_header_str

def test__section_header_to_str__no_args():
    slurm_configuration_list = [
        ['job_name', 'job_name'],
        ['qos', 'phillpot'],
        ['mail_type', 'END'],
        ['mail', 'eragasa@ufl.edu'],
        ['ntasks', 32],
        ['output', 'job.out'],
        ['error', 'job.err'],
        ['time','1:00:00'],
        ['memory', '4gb']
    ]

    slurm_configuration_dict = OrderedDict(
        [(k[0],k[1]) for k in slurm_configuration_list]
    )

    slurm_header_str = "\n".join([
        "#!/bin/bash",
        "#SBATCH --job-name=job_name",
        "#SBATCH --qos=phillpot",
        "#SBATCH --mail-type=END",
        "#SBATCH --mail-user=eragasa@ufl.edu",
        "#SBATCH --ntasks=32",
        "#SBATCH --distribution=cyclic:cyclic",
        "#SBATCH --mem=4gb",
        "#SBATCH --time=1:00:00",
        "#SBATCH --output=job.out",
        "#SBATCH --error=job.err"
    ]) + "\n"

    slurm_script = SlurmSubmissionScript()
    assert isinstance(slurm_script, SlurmSubmissionScript)
    assert isinstance(slurm_script.configuration, OrderedDict)

    slurm_script.configuration = slurm_configuration_dict
    assert isinstance(slurm_script.configuration, OrderedDict)
    slurm_header_section = slurm_script.section_header_section_to_str()
    print(slurm_header_section)
    print(slurm_header_str)
    assert slurm_header_section == slurm_header_str

if __name__ == "__main__":
    dev__section_header_to_str__no_args()
