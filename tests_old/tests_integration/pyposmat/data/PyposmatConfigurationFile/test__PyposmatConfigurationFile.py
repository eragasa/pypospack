import os
from collections import OrderedDict

# import from pypospack
import pypospack.utils
from pypospack.pyposmat.data import PyposmatConfigurationFile

# imports for the testing framework
import pytest

_pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()

_scenario_MgO_buckingham = OrderedDict()
_scenario_MgO_buckingham['data_input_dir'] = os.path.join(_pypospack_root_dir,'data','MgO_pareto_data')
_scenario_MgO_buckingham['data_output_dir'] = os.path.join(_pypospack_root_dir,'data','MgO_pareto_data')
_scenario_MgO_buckingham['configuration_fn'] = os.path.join(_scenario_MgO_buckingham['data_input_dir'],'pyposmat.config.in')
_scenario_MgO_buckingham['parameter_names'] = ['chrg_Mg', 'chrg_O', 'MgMg_A', 'MgMg_rho', 'MgMg_C', 'MgO_A', 'MgO_rho', 'MgO_C', 'OO_A', 'OO_rho', 'OO_C']
_scenario_MgO_buckingham['free_parameter_names'] = ['chrg_Mg', 'MgO_A', 'MgO_rho', 'OO_A', 'OO_rho', 'OO_C']

def test__init__no_args():
    o_config = PyposmatConfigurationFile()

def test__init__filename_args():
    _config_fn = _scenario_MgO_buckingham['configuration_fn']
    _parameter_names = _scenario_MgO_buckingham['parameter_names']
    _free_parameter_names = _scenario_MgO_buckingham['free_parameter_names']
    
    o_config = PyposmatConfigurationFile(filename=_config_fn)

    # check parameter_names, order doesn't matter
    # 1. check that the lists are of the same size
    # 2. check that if x in A, then x also in B
    
    # parameter check task 1
    assert len(o_config.parameter_names) == len(_parameter_names)

    # parameter check task 2
    for pn in o_config.parameter_names:
        assert pn in _parameter_names

    # check free parameter_names, order doesn't matter
    # 1. check that the list are of the same size
    # 2. check that if x in A, then x also in B
    assert len(o_config.free_parameter_names) == len(_free_parameter_names)
    for pn in o_config.free_parameter_names:
        assert pn in _free_parameter_names

def test__read_configuration_file():
    _config_fn = _scenario_MgO_buckingham['configuration_fn']
    _parameter_names = _scenario_MgO_buckingham['parameter_names']
    _free_parameter_names = _scenario_MgO_buckingham['free_parameter_names']
    o_config = PyposmatConfigurationFile()
    o_config.read(filename=_config_fn)

    # check parameter_names, order doesn't matter
    # 1. check that the lists are of the same size
    # 2. check that if x in A, then x also in B
    
    # parameter check task 1
    assert len(o_config.parameter_names) == len(_parameter_names)

    # parameter check task 2
    for pn in o_config.parameter_names:
        assert pn in _parameter_names

    # check free parameter_names, order doesn't matter
    # 1. check that the list are of the same size
    # 2. check that if x in A, then x also in B
    assert len(o_config.free_parameter_names) == len(_free_parameter_names)
    for pn in o_config.free_parameter_names:
        assert pn in _free_parameter_names

if __name__ == "__main__":
    _config_fn = _scenario_MgO_buckingham['configuration_fn']
    o_config = PyposmatConfigurationFile()
    o_config.read(filename=_config_fn)

    #parameter_names
    print('parameter_names='+str([k for k in o_config.configuration['sampling_dist']]))

    #free parameter names
    print('free_parameter_names='+str([k for k,v in o_config.configuration['sampling_dist'].items() if v[0] != 'equals']))

