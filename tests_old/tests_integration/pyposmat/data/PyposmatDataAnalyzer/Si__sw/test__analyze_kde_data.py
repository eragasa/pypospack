import pytest
import os
from collections import OrderedDict
import yaml
import pypospack.utils
from pypospack.io.filesystem import OrderedDictYAMLLoader
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatDataAnalyzer

def get_testing_set():

    testing_set = OrderedDict()
    testing_set['data_directory'] = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'data','Si__sw__data','pareto_optimization_1'
    )
    testing_set['config_fn'] = os.path.join(
            testing_set['data_directory'],
           'pyposmat.config.in'
    )
    testing_set['analysis_fn'] = 'pyposmat.analysis.out'

    assert os.path.isdir(testing_set['data_directory'])
    assert os.path.isfile(testing_set['config_fn'])

    return testing_set

def test__analyze_kde_data__w_filename():

    testing_set = get_testing_set()
    i_iteration = 1

    kde_fn = os.path.join(
            testing_set['data_directory'],
            'pyposmat.kde.{}.out'.format(i_iteration+1)
    )
        
    o = PyposmatDataAnalyzer()
    o.initialize_configuration(config_fn=testing_set['config_fn'])

    o.analyze_kde_data(
            i_iteration=i_iteration,
            filename=kde_fn
    )

    assert o.kde_data_fn == kde_fn
    assert isintance(o.kde_statistics,OrderedDict)

    descriptive_statistics_categories = ['qoi']

    for k in descriptive_statistics_categories:
        assert k in o.kde_statistics

def dev__analyze_kde_data():

    testing_set = get_testing_set()
    i_iteration = 1

    kde_fn = os.path.join(
            testing_set['data_directory'],
            'pyposmat.kde.{}.out'.format(i_iteration+1)
    )
    print('args:')
    print('\tkde_fn:{}'.format(kde_fn))
    print('\tos.path.isfile(kde_fn):{}'.format(os.path.isfile(kde_fn)))
        
    o = PyposmatDataAnalyzer()
    o.initialize_configuration(config_fn=testing_set['config_fn'])

    o.analyze_kde_data(
            i_iteration=i_iteration,
            filename=kde_fn
    )

    print("type(o.kde_statistics:{}".format(str(type(o.kde_statistics))))
    for k in o.kde_statistics:
        print(k)
   
    print(o.str__descriptive_statistics(o.kde_statistics))


if __name__ == "__main__":
    dev__analyze_kde_data()
