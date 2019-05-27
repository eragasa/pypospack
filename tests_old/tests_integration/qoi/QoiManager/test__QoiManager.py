import pytest

import os,shutil
from collections import OrderedDict


import pypospack.utils
from pypospack.qoi import QoiManager
from pypospack.pyposmat.data import PyposmatConfigurationFile
pypospack_root_directory = pypospack.utils.get_pypospack_root_directory()
config_fn = os.path.join(pypospack_root_directory,'data/MgO_pareto_data/pyposmat.config.in')
simulation_directories = [
        'MgO_NaCl.lmps_min_all','MgO_NaCl.lmps_elastic','MgO_NaCl_001s.lmps_min_pos',
        'MgO_NaCl_fr_a.lmps_min_pos','MgO_NaCl_fr_c.lmps_min_pos','MgO_NaCl_sch.lmps_min_pos']

def test__import_class():
    from pypospack.qoi import QoiManager

def get_qoi_database_from_PyposmatConfigurationFile(config_fn):
    config = PyposmatConfigurationFile()
    config.read(filename=config_fn)

    assert type(config.qois) is OrderedDict
    for qoi_id, qoi_info in config.qois.items():
        assert type(qoi_id) is str
        assert set(qoi_info.keys()) == set(['qoi_type','structures','target'])

    return config.qois

def get_TaskManager(config_fn);
    
def get_parameter_set():
    from pypospack.potentials.MgO import MgO_LewisCatlow
    parameters=MgO_LewisCatlow['parameters']
    return parameters

def test____init____fullauto_True__qoi_database_OrderedDict():
    qoi_database=get_qoi_database_from_PyposmatConfigurationFile(config_fn=config_fn)

    o = QoiManager(qoi_database=qoi_database,fullauto=True)

def test____init____fullauto_False__qoi_database_OrderedDict():
    qoi_database=get_qoi_database_from_PyposmatConfigurationFile(config_fn=config_fn)

    o = QoiManager(qoi_database=qoi_database,fullauto=False)

def dev____init__():
    qoi_database = get_qoi_database_from_PyposmatConfigurationFile(config_fn=config_fn)

    print('config.qois')
    print('\ttype:',type(qoi_database))
    for qoi_id,qoi_info in qoi_database.items():
        print('\t{}:'.format(qoi_id))
        for k,v in qoi_info.items():
            print('\t\t',k,':',v)
if __name__ == "__main__":
    dev____init__()
    


