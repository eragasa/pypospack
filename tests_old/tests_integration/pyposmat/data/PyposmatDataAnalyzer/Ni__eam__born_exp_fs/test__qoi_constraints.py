import pytest
from collections import OrderedDict
import os

import pypospack.utils
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data.data_analyzer import NewPyposmatDataAnalyzer

class PyposmatDataAnalyzer(NewPyposmatDataAnalyzer):pass

def get_testing_set():
    testing_set = OrderedDict()
    testing_set['i_iteration'] = 4
    testing_set['config_fn']= os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.config.in'
    )
    testing_set['results_fn'] = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.results.{}.out'.format(
                testing_set['i_iteration']
            )
    )
    return testing_set

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
    o.analyze_results(i_iteration=i_iteration)

def dev__get_descriptive_statistics():
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

def dev__get_weights_by_z_error_normalization():
    testing_set = get_testing_set()
    config_fn = testing_set['config_fn']
    results_data_fn = testing_set['results_fn']

    assert os.path.isfile(config_fn)
    assert os.path.isfile(results_data_fn)

    o = PyposmatDataAnalyzer(config_fn=config_fn,results_data_fn=results_data_fn)
    weights = o.get_weights_by_z_error_normalization()

    assert isinstance(weights,OrderedDict)
    for k,v in weights.items():
        assert k in o.qoi_names
        assert isinstance(v,float)

    print(o.str__cost_function_weights(weights=weights))

def dev_get_weights_by_qoi_target_normalization():
    testing_set = get_testing_set()
    config_fn = testing_set['config_fn']
    results_data_fn = testing_set['results_fn']

    assert os.path.isfile(config_fn)
    assert os.path.isfile(results_data_fn)

    o = PyposmatDataAnalyzer(config_fn=config_fn,results_data_fn=results_data_fn)
    weights = o.get_weights_by_qoi_target_normalization()

    assert isinstance(weights,OrderedDict)
    for k,v in weights.items():
        assert k in o.qoi_names
        assert isinstance(v,float)
   
    print(o.str__cost_function_weights(weights=weights))

def test__filter_by_pareto_set_membership():
    testing_set = get_testing_set()
    config_fn = testing_set['config_fn']
    results_data_fn = testing_set['results_fn']

    assert os.path.isfile(config_fn)
    assert os.path.isfile(results_data_fn)

    o = PyposmatDataAnalyzer(config_fn=config_fn,results_data_fn=results_data_fn)
    is_survive_idx,pareto_set_info = o.filter_by_pareto_membership()

    assert all([isinstance(v,int) for v in is_survive_idx])
    assert isinstance(is_survive_idx,set)
    assert isinstance(pareto_set_info,OrderedDict)

def dev__filter_by_pareto_set_membership():
    print(80*'-')
    print('{:^80}'.format('filter_by_pareto_membership'))
    
    testing_set = get_testing_set()
    config_fn = testing_set['config_fn']
    results_data_fn = testing_set['results_fn']

    assert os.path.isfile(config_fn)
    assert os.path.isfile(results_data_fn)

    o = PyposmatDataAnalyzer(config_fn=config_fn,results_data_fn=results_data_fn)
    is_survive_idx,pareto_set_info = o.filter_by_pareto_membership()

    print(is_survive_idx)
    assert all([isinstance(v,int) for v in is_survive_idx])
    assert isinstance(is_survive_idx,set)
    assert isinstance(pareto_set_info,OrderedDict)
    print('type(is_survive_idx):{}'.format(type(is_survive_idx)))
    print('type(pareto_set_info):{}'.format(type(pareto_set_info)))

def dev__filter_by_cost_function():
    print(80*'-')
    print('{:^80}'.format('filter_by_cost_functio'))
    
    testing_set = get_testing_set()
    config_fn = testing_set['config_fn']
    results_data_fn = testing_set['results_fn']

    assert os.path.isfile(config_fn)
    assert os.path.isfile(results_data_fn)

    o = PyposmatDataAnalyzer(config_fn=config_fn,results_data_fn=results_data_fn)
    weights = o.get_weights_by_z_error_normalization()
    pct_to_keep = 0.95
    is_survive_idx,filter_info = o.filter_by_cost_function(
            loss_function_type='abs_error',
            cost_function_type='weighted_sum',
            weights=weights,
            pct_to_keep=pct_to_keep
    )
    print('type(is_survive_idx):{}'.format(type(is_survive_idx)))
    print('type(cost_function_info:{}'.format(type(filter_info)))

    assert isinstance(is_survive_idx,set)
    assert isinstance(filter_info,OrderedDict)

    print('<---- START: str__filter_by_cost_function_filter_info')
    print(o.str__filter_by_cost_function_filter_info(filter_info=filter_info))
    print('<---- END: str__filter_by_cost_function_filter_info')

def test__filter_by_qoi_constraints():
    testing_set = get_testing_set()
    config_fn = testing_set['config_fn']
    results_data_fn = testing_set['results_fn']

    assert os.path.isfile(config_fn)
    assert os.path.isfile(results_data_fn)

    o = PyposmatDataAnalyzer(config_fn=config_fn,results_data_fn=results_data_fn)
    is_survive_idx,qoi_constraint_info = o.filter_by_qoi_constraints()

    assert isinstance(is_survive_idx,set)
    assert isinstance(qoi_constraint_info,OrderedDict)

def dev__filter_by_qoi_constraints():
    print(80*'-')
    print('{:^80}'.format('filter_by_qoi_constraints'))
    testing_set = get_testing_set()
    config_fn = testing_set['config_fn']
    results_data_fn = testing_set['results_fn']

    assert os.path.isfile(config_fn)
    assert os.path.isfile(results_data_fn)
    o = PyposmatDataAnalyzer(config_fn=config_fn,results_data_fn=results_data_fn)

    #for k,v in o.qoi_constraints.items():
    #    print(k,v)
    #for k,v in o.configuration.qoi_constraints.items():
    #    print(k,v)

    is_survive_idx,qoi_constraint_info = o.filter_by_qoi_constraints()
    
    print('type(is_survive_idx):{}'.format(type(is_survive_idx)))
    assert isinstance(is_survive_idx,set)
    print('type(qoi_constraint_info:{}'.format(type(qoi_constraint_info)))
    assert isinstance(qoi_constraint_info,OrderedDict)
    print(o.str__qoi_constraint_info(qoi_constraint_info=qoi_constraint_info))
if __name__ == "__main__":
    #dev__filter_by_cost_function()
    #dev__property__n_potentials_start()
    #dev__get_descriptive_statistics()
    #dev__get_weights_by_z_error_normalization()
    dev_get_weights_by_qoi_target_normalization()
    #dev__analyze_results()
    #dev__filter_by_pareto_set_membership()
    dev__filter_by_qoi_constraints()
