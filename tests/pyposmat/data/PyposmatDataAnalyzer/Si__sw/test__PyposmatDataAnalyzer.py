import pytest
import os
import pandas as pd
from collections import OrderedDict
import pypospack.utils
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data.data_analyzer import PyposmatDataAnalyzer

def get_testing_set():
    testing_set = OrderedDict()
    testing_set['config_fn'] = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'data','Si__sw__data','pareto_optimization_1','pyposmat.config.in')
    testing_set['results_data_fn'] = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'data','Si__sw__data','pareto_optimization_1','pyposmat.results.9.out')

    assert os.path.isfile(testing_set['config_fn'])
    assert os.path.isfile(testing_set['results_data_fn'])

    return testing_set

def test____init____no_arguments():
    o = PyposmatDataAnalyzer()

    assert type(o) is PyposmatDataAnalyzer
    assert o.configuration is None
    assert o.results_data is None
    assert o.results_data_fn is None
    assert o.config_fn is None

def test____init____using_path_args():
    testing_set = get_testing_set()

    o = PyposmatDataAnalyzer(
            config_fn=testing_set['config_fn'],
            results_data_fn=testing_set['results_data_fn']
    )

    config = PyposmatConfigurationFile()
    config.read(filename=testing_set['config_fn'])

    assert isinstance(o,PyposmatDataAnalyzer)

    assert o.config_fn == testing_set['config_fn']
    assert isinstance(o.configuration,PyposmatConfigurationFile)
    assert isinstance(o.parameter_names,list)
    assert set(o.parameter_names) == set(config.parameter_names)
    assert isinstance(o.error_names, list)
    assert set(o.error_names) == set(config.error_names)
    assert isinstance(o.qoi_names, list)
    assert set(o.qoi_names) == set(config.qoi_names)
    
    assert o.results_data_fn == testing_set['results_data_fn']
    assert isinstance(o.results_data,PyposmatDataFile)
    assert isinstance(o.results_df,pd.DataFrame)

def test__initialize_configuration__no_args():
    testing_set = get_testing_set()
    config_fn = testing_set['config_fn']
    assert isinstance(testing_set['config_fn'],str)

    o = PyposmatDataAnalyzer()
    o.config_fn = config_fn
    o.initialize_configuration()

    assert o.config_fn == testing_set['config_fn']
    assert type(o.configuration) is PyposmatConfigurationFile

def test__initialize_configuration__no_args__no_attributes():
    o = PyposmatDataAnalyzer()

    with pytest.raises(TypeError) as e:
        o.initialize_configuration()

def test__initialize_configuration__with_path():
    testing_set = get_testing_set()
    assert isinstance(testing_set['config_fn'],str)

    o = PyposmatDataAnalyzer()
    o.initialize_configuration(config_fn=testing_set['config_fn'])

def dev__initialize_configuration__with_path():
    print(80*'-')
    print('{:^80}'.format('initialize_configuration__with_path'))
    testing_set = get_testing_set()
    o = PyposmatDataAnalyzer()
    o.initialize_configuration(config_fn=testing_set['config_fn'])

def test__initialize_configuration__with_object():
    testing_set = get_testing_set()

    o_config = PyposmatConfigurationFile()
    o_config.read(filename=testing_set['config_fn'])

    o = PyposmatDataAnalyzer()
    o.initialize_configuration(o_config=o_config)

def test__initialize_configuration__with_object_and_path():
    testing_set = get_testing_set()

    o_config = PyposmatConfigurationFile()
    o_config.read(filename=testing_set['config_fn'])

    o = PyposmatDataAnalyzer()
    with pytest.raises(TypeError) as e:
        o.initialize_configuration(config_fn=testing_set['config_fn'],
                                  o_config=o_config)

def test__property__n_potentials_start():
    testing_set = get_testing_set()
    config_fn = testing_set['config_fn']
    results_data_fn = testing_set['results_fn']

    o = PyposmatDataAnalyzer(config_fn=config_fn,results_data_fn=results_data_fn)

    assert isinstance(o.n_potentials_start,int)

def dev__property__n_potentials_start():
    print(80*'-')
    print('{:^80}'.format('property -> n_potentials_start'))
    testing_set = get_testing_set()
    config_fn = testing_set['config_fn']
    results_data_fn = testing_set['results_fn']

    o = PyposmatDataAnalyzer(config_fn=config_fn,results_data_fn=results_data_fn)

    print("type(o.n_potentials_start)):{}".format(type(o.n_potentials_start)))
    print("o.n_potentials_start:{}".format(o.n_potentials_start))

if __name__ == "__main__":
    #dev__initialize_configuration__with_path()
    #dev__read_configuration_file()
    dev__qoi_constraints()
