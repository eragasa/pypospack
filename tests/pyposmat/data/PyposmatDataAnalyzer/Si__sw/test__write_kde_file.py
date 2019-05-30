import pytest
import os
from collections import OrderedDict
import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataAnalyzer

def get_testing_set():
    testing_set = OrderedDict()
    testing_set['i_iteration'] = 1
    testing_set['config_fn'] = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'data','Si__sw__data','pareto_optimization_1','pyposmat.config.in'
    )
    testing_set['results_data_fn'] = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'data','Si__sw__data','pareto_optimization_1',
            'pyposmat.results.{}.out'.format(testing_set['i_iteration']-1)
    )
    
    assert os.path.isfile(testing_set['config_fn'])
    assert os.path.isfile(testing_set['results_data_fn'])

    return testing_set

def dev__write_kde_file():
    kde_data_fn = 'pyposmat.kde.test.out'

    try:
        os.remove(kde_data_fn)
    except FileNotFoundError as e:
        pass

    testing_set = get_testing_set()

    o = PyposmatDataAnalyzer()
    o.initialize_configuration(config_fn=testing_set['config_fn'])
    o.initialize_results_data(results_data_fn=testing_set['results_data_fn'])
    is_survive_idx,filter_info = o.analyze_results(i_iteration=testing_set['i_iteration'])

    o.write_kde_file(filename=kde_data_fn)

    print(o.kde_data.df['sim_id'])
if __name__ == "__main__":
    dev__write_kde_file()
