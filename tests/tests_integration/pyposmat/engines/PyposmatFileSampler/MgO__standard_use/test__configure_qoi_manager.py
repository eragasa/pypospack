import os
from collections import OrderedDict
import pytest

import pypospack.utils
from pypospack.qoi import QoiManager
from pypospack.pyposmat.engines import PyposmatFileSampler
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

test_case_MgO = OrderedDict()
test_case_MgO['config_fn'] = os.path.join(
        pypospack.utils.get_pypospack_root_directory(),
        'data/MgO_pareto_data/pyposmat.config.in'
        )
test_case_MgO['data_in_fn'] = os.path.join(
        pypospack.utils.get_pypospack_root_directory(),
        'data/MgO_pareto_data/culled_004.out'
        )
test_case_MgO['qoi_names'] = [
        'MgO_NaCl.a0',
        'MgO_NaCl.c11', 'MgO_NaCl.c12', 'MgO_NaCl.c44', 
        'MgO_NaCl.fr_a', 'MgO_NaCl.fr_c', 'MgO_NaCl.sch', 
        'MgO_NaCl.B', 'MgO_NaCl.G', 'MgO_NaCl.001s']
def show_test_case(test_case):
    """show the results of the test case

    Args:
        test_case(OrderedDict): a test case fixture

    """

    assert isinstance(test_case,OrderedDict)

    s = 'config_fn={}'.format(test_case['config_fn'])
    print(s)
    s = 'data_in_fn={}'.format(test_case['data_in_fn'])
    print(s)

def test__configure_qoi_manager():
    test_case = test_case_MgO

    o = PyposmatFileSampler(config_fn=test_case['config_fn'],
                            data_in_fn=test_case['data_in_fn'],
                            fullauto=False)
    o.read_configuration_file(filename=test_case['config_fn'])
    o.configure_qoi_manager() 

    assert isinstance(o.configuration,PyposmatConfigurationFile)
    assert isinstance(o.qoi_manager,QoiManager)

def dev__configure_qoi_manager():
    config_fn = test_case_MgO['config_fn']
    data_in_fn = test_case_MgO['data_in_fn']

    o = PyposmatFileSampler(config_fn=config_fn,
                            data_in_fn=data_in_fn) 
    qois = [k for k in o.configuration.qois]
    print('qoi_names:',qois)
    o.configure_qoi_manager()
    print('o.configuration.qoi_names:',o.configuration.qoi_names)
    s = str(type(o.qoi_manager))
    print("type(o.qoi_manager):{}".format(s))
    
    s = str(type((o.qoi_manager.tasks)))
    print('type(o.qoi_manager.tasks:{}'.format(s))

if __name__ == "__main__":
    show_test_case(test_case_MgO)
    dev__configure_qoi_manager()

