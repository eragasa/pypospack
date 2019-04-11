import pytest
import os
from collections import OrderedDict

import pypospack.utils
from pypospack.pyposmat.data import PyposmatBadParametersFile

def get_testing_set():
    pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
    testing_set = OrderedDict()
    testing_set['data_directory'] = os.path.join(
            pypospack_root_dir,'examples','Si__sw',
            'dev__pareto_optimization','data')
    testing_set['config_fn'] = os.path.join(
            testing_set['data_directory'],'pyposmat.config.in')
    testing_set['badparameters_out_fn'] = os.path.join('pyposmat.badparameters.out')
    testing_set['badparameters_in_fn'] = os.path.join('pyposmat.badparameters.in')

    return testing_set

def test____init____wo_filename():
    o = PyposmatBadParametersFile()
    assert o.filename == 'pyposmat.badparameters.out'
    assert o.configuration is None

def test____init____w_filename():
    test_badparameters_fn = 'test_filename.out'

    if os.path.isfile(test_badparameters_fn):
        os.remove(test_badparameters_fn)

    o = PyposmatBadParametersFile(filename=test_badparameters_fn)
    assert o.filename == test_badparameters_fn
    assert not os.path.isfile(test_badparameters_fn)
    assert o.configuration is None

def test____init____w_filename_config_fn():
    testing_set = get_testing_set()

    assert os.path.isfile(testing_set['badparameters_in_fn'])
    assert os.path.isfile(testing_set['config_fn'])

    if os.path.isfile(testing_set['badparameters_out_fn']):
        m = "removing the badparameter_out_file:{}".format(
                testing_set['badparameters_out_fn'])
        print(m)
        os.remove(testing_set['badparameters_out_fn'])

    o = PyposmatBadParametersFile(
            filename=testing_set['badparameters_out_fn'],
            config_fn=testing_set['config_fn'])

    assert o.filename == testing_set['badparameters_out_fn']

    from pypospack.pyposmat.data import PyposmatConfigurationFile
    o_config = PyposmatConfigurationFile()
    o_config.read(filename=testing_set['config_fn'])
    assert isinstance(o.configuration,PyposmatConfigurationFile)
    assert isinstance(o.parameter_names,list)
    assert set(o.parameter_names) == set(o_config.parameter_names)
