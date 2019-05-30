"""
This module tests the PyposmatEngine for a buckingham simulation of MgO using the Lewis and Catlow potential
"""
import pytest
import os,shutil
from collections import OrderedDict
import pypospack.utils
from pypospack.pyposmat.engines import PyposmatEngine
from pypospack.pyposmat.data import PyposmatDataFile


pypospack_root_directory = pypospack.utils.get_pypospack_root_directory()
config_fn = os.path.join(pypospack_root_directory,'data/MgO_pareto_data/pyposmat.config.in')
simulation_directories = [
        'MgO_NaCl.lmps_min_all','MgO_NaCl.lmps_elastic','MgO_NaCl_001s.lmps_min_pos',
        'MgO_NaCl_fr_a.lmps_min_pos','MgO_NaCl_fr_c.lmps_min_pos','MgO_NaCl_sch.lmps_min_pos']

from pypospack.potentials.MgO import MgO_LewisCatlow
parameters=MgO_LewisCatlow['parameters']

def cleanup_simulation_directories():
    for d in simulation_directories:
        try:
           shutil.rmtree(d)
        except FileNotFoundError as e:
            pass

def test__init__():
    engine = PyposmatEngine(filename_in=config_fn,filename_out='pyposmat.config.out')

def test__create_base_directories():
    engine = PyposmatEngine(filename_in=config_fn,filename_out='pyposmat.config.out')
    engine.create_base_directories()

    assert all([os.path.isdir(d) for d in simulation_directories])

    cleanup_simulation_directories()

def test__read_configuration_file():
    engine = PyposmatEngine(filename_in=config_fn,filename_out='pyposmat.config.out')
    engine.create_base_directories()
    engine.read_configuration_file()

    from pypospack.pyposmat.data import PyposmatConfigurationFile
    assert type(engine.configuration) is PyposmatConfigurationFile

    cleanup_simulation_directories()

def test__configure_qoi_manager():
    engine = PyposmatEngine(filename_in=config_fn,filename_out='pyposmat.config.out')
    engine.create_base_directories()
    engine.read_configuration_file()
    engine.configure_qoi_manager()

    from pypospack.qoi import QoiManager
    assert type(engine.qoi_manager) is QoiManager

    cleanup_simulation_directories()

def test__configure_task_manager():
    engine = PyposmatEngine(filename_in=config_fn,filename_out='pyposmat.config.out')
    engine.create_base_directories()
    engine.read_configuration_file()
    engine.configure_qoi_manager()
    engine.configure_task_manager()

    from pypospack.task import TaskManager
    assert type(engine.task_manager) is TaskManager

    cleanup_simulation_directories()
def test__configure():
    engine = PyposmatEngine(filename_in=config_fn,filename_out='pyposmat.config.out')
    engine.configure()
    cleanup_simulation_directories()

def test__evaluate_parameter_set():
    engine = PyposmatEngine(filename_in=config_fn,filename_out='pyposmat.config.out')
    engine.configure()
    results = engine.evaluate_parameter_set(parameters=parameters)
    cleanup_simulation_directories()
