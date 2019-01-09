import pytest
import os,shutil
from collections import OrderedDict
import pypospack.utils
from pypospack.pyposmat.engines import PyposmatEngine
from pypospack.pyposmat.data import PyposmatDataFile

_pypospack_root_directory = pypospack.utils.get_pypospack_root_directory()
_pyposmat_config_fn = os.path.join(_pypospack_root_directory,'data','MgO_pareto_data','pyposmat.config.in')

_simulation_directories = [
        'MgO_NaCl.lmps_min_all',
        'MgO_NaCl.lmps_elastic',
        'MgO_NaCl_001s.lmps_min_pos',
        'MgO_NaCl_fr_a.lmps_min_pos',
        'MgO_NaCl_fr_c.lmps_min_pos',
        'MgO_NaCl_sch.lmps_min_pos']

def cleanup_simulation_directories():
    for d in _simulation_directories:
        shutil.rmtree(d)

def test____init____w_filenames():
    engine = PyposmatEngine(
            filename_in = _pyposmat_config_fn,
            filename_out = 'pypospack.config.out')
    engine.configure()

    for d in _simulation_directories:
        assert os.path.isdir(d)

    cleanup_simulation_directories()

def test____evaluate_parameter_set():
    from pypospack.potentials.MgO import MgO_LewisCatlow
    _parameters=MgO_LewisCatlow['parameters']

    engine = PyposmatEngine(
            filename_in=_pyposmat_config_fn,
            filename_out='pypospack.config.out')
    engine.configure()

    _results=engine.evaluate_parameter_set(parameters=_parameters)
    _results_parameters = _results['parameters']
    for k,v in _results_parameters.items():
        assert v is not None

    _results_qois = _results['qois']
    for k,v in _results_qois.items():
        assert v is not None

    cleanup_simulation_directories()
