import pytest

import os
import pypospack.utils
from pypospack.pyposmat.engines import PyposmatEngine
from pypospack.pyposmat.data import PyposmatDataFile

pypospack_root_directory = pypospack.utils.get_pypospack_root_directory()

config_fn = os.path.join(pypospack_root_directory,'data/MgO_pareto_data/pyposmat.config.in')

simulation_directories = [
        'MgO_NaCl.lmps_min_all',
        'MgO_NaCl.lmps_elastic',
        'MgO_NaCl_001s.lmps_min_pos',
        'MgO_NaCl_fr_a.lmps_min_pos',
        'MgO_NaCl_fr_c.lmps_min_pos',
        'MgO_NaCl_sch.lmps_min_pos']

from pypospack.potentials.MgO import MgO_LewisCatlow
parameters=MgO_LewisCatlow['parameters']
def test__configure_task_manager():
    engine = PyposmatEngine(filename_in=config_fn,filename_out='pyposmat.config.out')
    
    # these need to be done first
    engine.create_base_directories()
    engine.read_configuration_file()
    engine.configure_qoi_manager()

    # this line is what is tested
    engine.configure_task_manager()
    
    from pypospack.task import TaskManager
    assert type(engine.task_manager) is TaskManager
