import pytest
from collections import OrderedDict
import os

import pypospack.utils
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatDataAnalyzer
from pypospack.pyposmat.data.data_analyzer import PyposmatDataAnalyzer

def get_testing_set():

    testing_set = OrderedDict()
    testing_set['i_iteration'] = 1
    testing_set['config_fn']= os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'examples/Si__sw/alt_pareto_optimization/data',
            'pyposmat.config.in'
    )
    testing_set['results_fn'] = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'examples/Si__sw/alt_pareto_optimization/data',
            'pyposmat.results.{}.out'.format(testing_set['i_iteration'])
    )
    return testing_set

def test__analyze_results():

    testing_set = get_testing_set()
    i_iteration = testing_set['i_iteration']
    config_fn = testing_set['config_fn']
    results_data_fn = testing_set['results_fn']

    assert os.path.isfile(config_fn)
    assert os.path.isfile(results_data_fn)

    o = PyposmatDataAnalyzer(config_fn=config_fn,results_data_fn=results_data_fn)
    is_survive_idx,filter_info = o.analyze_results(i_iteration=i_iteration)

    assert isinstance(is_survive_idx,set)
    assert isinstance(filter_info,OrderedDict)

def dev__analyze_results():
    print(80*'-')
    print('{:^80}'.format('method -> analyze_results'))

    testing_set = get_testing_set()
    i_iteration = testing_set['i_iteration']
    config_fn = testing_set['config_fn']
    results_data_fn = testing_set['results_fn']

    assert os.path.isfile(config_fn)
    assert os.path.isfile(results_data_fn)

    o = PyposmatDataAnalyzer(config_fn=config_fn,results_data_fn=results_data_fn)
    is_survive_idx,filter_info = o.analyze_results(i_iteration=i_iteration)

    print('type(is_survive_idx):{}'.format(str(type(is_survive_idx))))
    print('type(filter_info):{}'.format(str(type(filter_info))))

if __name__ == "__main__":
    dev__analyze_results()
