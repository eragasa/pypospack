import pytest
from collections import OrderedDict
from pypospack.io.slurm import SlurmSubmissionScript

slurm_vasp_modules = [
    'intel/2016.0.109',
    'impi/5.1.1'
]

slurm_load_module_str = "module load intel/2016.0.109\n"
slurm_load_module_str += "module load impi/5.1.1\n"

def dev__section_load_modules_to_str__module_arg():
    slurm_script = SlurmSubmissionScript()
    assert isinstance(slurm_script, SlurmSubmissionScript)
    assert isinstance(slurm_script.configuration, OrderedDict)

    str_load_modules = slurm_script.section_load_modules_to_str(modules=slurm_vasp_modules)
    assert slurm_script.configuration['modules'] == slurm_vasp_modules

    assert str_load_modules == slurm_load_module_str
    print(str_load_modules)

def dev__section_load_modules_to_str__no_args():
    slurm_script = SlurmSubmissionScript()
    assert isinstance(slurm_script, SlurmSubmissionScript)
    assert isinstance(slurm_script.configuration, OrderedDict)

    slurm_script.configuration['modules'] = slurm_vasp_modules
    str_load_modules = slurm_script.section_load_modules_to_str()
    assert str_load_modules == slurm_load_module_str
    print(str_load_modules)

def test__section_load_modules_to_str__module_arg():
    slurm_script = SlurmSubmissionScript()
    assert isinstance(slurm_script, SlurmSubmissionScript)
    assert isinstance(slurm_script.configuration, OrderedDict)

    str_load_modules = slurm_script.section_load_modules_to_str(modules=slurm_vasp_modules)
    assert slurm_script.configuration['modules'] == slurm_vasp_modules

    assert str_load_modules == slurm_load_module_str

def test__section_load_modules_to_str__no_args():
    slurm_script = SlurmSubmissionScript()
    assert isinstance(slurm_script, SlurmSubmissionScript)
    assert isinstance(slurm_script.configuration, OrderedDict)

    slurm_script.configuration['modules'] = slurm_vasp_modules
    str_load_modules = slurm_script.section_load_modules_to_str()
    assert str_load_modules == slurm_load_module_str

if __name__ == "__main__":
    dev__section_load_modules_to_str__no_args()
    dev__section_load_modules_to_str__module_arg()
