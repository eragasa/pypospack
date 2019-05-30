import os
from collections import OrderedDict
import pytest

import pypospack.utils
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
test_case_MgO['i_iteration'] = 0

def show_test_case(test_case):
    """show the results of the test case

    Args:
        test_case(OrderedDict): a test case fixture

    """
    assert isinstance(test_case,OrderedDict)

    s = 'config_fn={}'.format(test_case['config_fn'])
    print(s)

def test__run_simulations():
    test_case = test_case_MgO

    assert os.path.isfile(test_case['config_fn'])
    assert os.path.isfile(test_case['data_in_fn'])

    o = PyposmatFileSampler(config_fn=test_case['config_fn'],
                            data_in_fn=test_case['data_in_fn'])
    o.run_simulations(i_iteration=test_case['i_iteration'],
                      filename=test_case['data_in_fn'])


def dev__run_simulations():
    test_case = test_case_MgO

    assert os.path.isfile(config_fn)
    assert os.path.isfile(data_in_fn)

    o = PyposmatFileSampler(config_fn=test_case['config_fn'],
                            data_in_fn=test_case['data_in_fn'])
    o.run_simulations(i_iteration=test_case['i_iteration'])

if __name__ == "__main__":
    show_test_case(test_case=test_case_MgO)
    dev__configure_qoi_manager()
    dev__run_simulations()

