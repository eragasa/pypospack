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
    testing_set['kde_fn'] = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'examples/Si__sw/alt_pareto_optimization/data',
            'pyposmat.kde.{}.out'.format(testing_set['i_iteration'])
    )
    return testing_set

def test__get_descriptive_statistics():
    testing_set = get_testing_set()
    config_fn = testing_set['config_fn']
    results_data_fn = testing_set['results_fn']

    assert os.path.isfile(config_fn)
    assert os.path.isfile(results_data_fn)

    o = PyposmatDataAnalyzer(config_fn=config_fn,results_data_fn=results_data_fn)
    descriptive_statistics = o.get_descriptive_statistics()

    assert isinstance(descriptive_statistics,OrderedDict)
    assert 'qois' in descriptive_statistics
    for qn in o.qoi_names:
        assert qn in descriptive_statistics['qois']

    print(o.str__descriptive_statistics(descriptive_statistics=descriptive_statistics))

def dev__get_descriptive_statistics():
    print(80*'-')
    print('{:^80}'.format('method -> get_descriptive_statistics'))
    
    testing_set = get_testing_set()
    config_fn = testing_set['config_fn']
    results_data_fn = testing_set['results_fn']

    assert os.path.isfile(config_fn)
    assert os.path.isfile(results_data_fn)

    o = PyposmatDataAnalyzer(config_fn=config_fn,results_data_fn=results_data_fn)
    descriptive_statistics = o.get_descriptive_statistics()

    assert isinstance(descriptive_statistics,OrderedDict)
    assert 'qois' in descriptive_statistics
    for qn in o.qoi_names:
        assert qn in descriptive_statistics['qois']

    print(o.str__descriptive_statistics(descriptive_statistics=descriptive_statistics))
    print(o.results_data.df.shape)

def dev__get_descriptive_statistics__from_kde_file():
    print(80*'-')
    print('{:^80}'.format('method -> get_descriptive_statistics__from_kde_file'))

    testing_set = get_testing_set()
    config_fn = testing_set['config_fn']
    results_data_fn = testing_set['results_fn']
    kde_data_fn = testing_set['kde_fn']

    assert os.path.isfile(config_fn)
    assert os.path.isfile(results_data_fn)
    assert os.path.isfile(kde_data_fn)

    o = PyposmatDataAnalyzer(config_fn=config_fn,results_data_fn=results_data_fn)

    kde_data = PyposmatDataFile()
    kde_data.read(filename=kde_data_fn)

    descriptive_statistics = o.get_descriptive_statistics(df=kde_data.df)

    print(o.str__descriptive_statistics(descriptive_statistics=descriptive_statistics))
    print(kde_data.df.shape)

if __name__ == "__main__":
    dev__get_descriptive_statistics()
    dev__get_descriptive_statistics__from_kde_file()
