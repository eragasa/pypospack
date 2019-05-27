import os,shutil
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

def cleanup_simulation_directories():
    for d in simulation_directories:
        shutil.rmtree(d)

def test__initialize_w_filenames():
    engine = PyposmatEngine(filename_in=config_fn,filename_out='pyposmat.config.out')
    engine.configure()

    for d in simulation_directories:
        assert os.path.isdir(d)

    cleanup_simulation_directories()

