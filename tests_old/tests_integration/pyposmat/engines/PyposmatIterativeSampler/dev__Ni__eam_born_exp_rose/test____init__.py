import pytest
import os
import pypospack.utils
from pypospack.pyposmat.engines import PyposmatIterativeSampler
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile

_pypospack_root = pypospack.utils.get_pypospack_root_directory()
_simulation_root = os.path.join(
    _pypospack_root,
    'examples','Ni__eam__born_exp_rose','01_preconditioning_3.5NN')
_pyposmat_configuration_fn = os.path.join(
            _simulation_root,'data','pyposmat.config.in')

def test____init__():
    engine = PyposmatIterativeSampler(
            configuration_filename=_pyposmat_configuration_fn,
            log_fn='pyposmat.log')

def test__read_configuration_files():
    engine = PyposmatIterativeSampler(
            configuration_filename=_pyposmat_configuration_fn,
            log_fn='pyposmat.log')
    engine.read_configuration_file(filename=_pyposmat_configuration_fn)
    assert type(engine.configuration) is PyposmatConfigurationFile

def test__read_configuration_files__no_args():
    engine = PyposmatIterativeSampler(
            configuration_filename=_pyposmat_configuration_fn,
            log_fn='pyposmat.log')
    engine.read_configuration_file()
    assert type(engine.configuration) is PyposmatConfigurationFile

