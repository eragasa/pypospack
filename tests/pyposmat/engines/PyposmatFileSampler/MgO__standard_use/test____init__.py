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

def show_test_case(test_case):
    """show the results of the test case

    Args:
        test_case(OrderedDict): a test case fixture

    """
    assert isinstance(test_case,OrderedDict)

    s = 'config_fn={}'.format(test_case['config_fn'])
    print(s)

def test____init__():
    config_fn = test_case_MgO['config_fn']
    data_in_fn = test_case_MgO['data_in_fn']
    
    assert os.path.isfile(config_fn)
    assert os.path.isfile(data_in_fn)

    o = PyposmatFileSampler(config_fn=config_fn,
                            data_in_fn=data_in_fn)

    assert isinstance(o.configuration, PyposmatConfigurationFile)
    assert isinstance(o.datafile_in, PyposmatDataFile)
    assert isinstance(o.datafile_out, PyposmatDataFile)

def dev____init__():
    config_fn = test_case_MgO['config_fn']
    data_in_fn = test_case_MgO['data_in_fn']

    assert os.path.isfile(config_fn)
    assert os.path.isfile(data_in_fn)

    o = PyposmatFileSampler(config_fn=config_fn,
                            data_in_fn=data_in_fn)

    s = []
    s.append("type(o.configuration)={}".format(type(o.configuration)))
    s.append("type(o.datafile_in)={}".format(type(o.datafile_in)))
    s.append("type(o.datafile_out={}".format(type(o.datafile_out)))
    print("\n".join(s))
if __name__ == "__main__":
    show_test_case(test_case=test_case_MgO)
    dev____init__()

