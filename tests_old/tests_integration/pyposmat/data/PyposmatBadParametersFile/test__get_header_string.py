import pytest
import os
from collections import OrderedDict
import pypospack.utils 
from pypospack.pyposmat.data import PyposmatBadParametersFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

def get_testing_set():
    pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
    testing_set = OrderedDict()
    testing_set['config_fn'] = os.path.join(
        pypospack.utils.get_pypospack_root_directory(),
        'examples','Si__sw','dev__pareto_optimization',
        'data','pyposmat.config.in')
    testing_set['badparameters_out_fn'] = 'pyposmat.badparameters.out'
    return testing_set

def dev__get_header_string():
    testing_set = get_testing_set()

    print(os.path.isfile(testing_set['config_fn']),testing_set['config_fn'])

    o = PyposmatBadParametersFile(
            filename = testing_set['badparameters_out_fn'],
            config_fn = testing_set['config_fn'])
    print(o.get_header_string())

def test__get_header_string():
    testing_set = get_testing_set()
    
    o_config = PyposmatConfigurationFile()
    o_config.read(filename=testing_set['config_fn'])
    assert os.path.isfile(testing_set['config_fn'])

    f = PyposmatBadParametersFile(
            filename = testing_set['badparameters_out_fn'],
            config_fn=testing_set['config_fn'])
    s = f.get_header_string()

    assert type(s) is str

    header_line_1 = ['sim_id'] \
            + o_config.parameter_names \
            + ['reason']
    header_line_2 = ['sim_id'] \
            + len(o_config.parameter_names)*['param'] \
            + ['reason']
    s_test = "{}\n".format(",".join(header_line_1))
    s_test += "{}\n".format(",".join(header_line_2))

    assert s_test == s

def test__write_header_string__no_args():
    testing_set = get_testing_set()

    o = PyposmatBadParametersFile(
            filename = testing_set['badparameters_out_fn'],
            config_fn=testing_set['config_fn'])
    o.write_header_section()

    assert os.path.isfile(testing_set['badparameters_out_fn'])

    if os.path.isfile(testing_set['badparameters_out_fn']):
        os.remove(testing_set['badparameters_out_fn'])

def test__write_header_string__w_filename():
    testing_set = get_testing_set()
    badparameters_out_fn = 'pyposmat.badparameters.test'

    o = PyposmatBadParametersFile(
            filename = testing_set['badparameters_out_fn'],
            config_fn=testing_set['config_fn'])
    o.write_header_section(filename=badparameters_out_fn)

    assert os.path.isfile(badparameters_out_fn)

    if os.path.isfile(badparameters_out_fn):
        os.remove(badparameters_out_fn)

if __name__ == "__main__":
    dev__get_header_string()
