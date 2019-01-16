import os
import pandas as pd
import numpy as np
from collections import OrderedDict

import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatDataAnalyzer
from pypospack.pyposmat.data import PyposmatConfigurationFile

class PyposmatUnknownPotentialScoringMetric(Exception):pass

def get_best_parameterization(config_fn,data_fn,metric_name='d_metric',o_config=None,o_data=None):
    _analyzer = PyposmatDataAnalyzer()
    _analyzer.read_configuration_file(filename=config_fn)
    _analyzer.read_data_file(filename=data_fn)

    # calculate the scoring metric
    if metric_name is 'd_metric':
        _df = _analyzer.calculate_d_metric(df=_analyzer.datafile.df)
    else:
        s = "The metric name {} is unsupported"
        s = s.format(metric_name)
        raise PyposmatUnsupportedPotentialScoringMetric(s)

    _data = PyposmatDataFile()
    _data.read(filename=data_fn)
    _data.df = _df
    _data.subselect_by_score(score_name='d_metric',n=1)

    _free_parameter_names = _analyzer.configuration.free_parameter_names
    
    _parameter_best_dict = OrderedDict()
    for pn in _free_parameter_names:
        _parameter_best_dict[pn] = _data.sub_parameter_df.iloc[0][pn]

    return _parameter_best_dict

def get_parameter_variance(
        config_fn,data_fn,
        metric_name='d_metric',
        n=100,
        o_config=None,
        o_data=None):
    """
    Args:
        config_fn (str):
        data_fn (str):
        metric_name (str):  (default:d_metric)
        n (int): the number of best metric values
        o_config (pypospack.config.data.PyposmatConfigurationFile)
        o_data (pypospack.config.data.PyposmatDataFile)
    Returns:
        collections.OrderedDict
    Raises:
        PyposmatUnknownPotentialScoringMetric
    """

    _analyzer = PyposmatDataAnalyzer()
    _analyzer.read_configuration_file(filename=config_fn)
    _analyzer.read_data_file(filename=data_fn)
   
    # calculate the scoring metric
    if metric_name is 'd_metric':
        _df = _analyzer.calculate_d_metric(df=_analyzer.datafile.df)
    else:
        s = "The metric name {} is unsupported"
        s = s.format(metric_name)
        raise PyposmatUnsupportedPotentialScoringMetric(s)

    _data = PyposmatDataFile()
    _data.read(filename=data_fn)
    _data.df = _df
    _data.subselect_by_score(score_name='d_metric',n=n)

    _param_std_df = _data.sub_parameter_df.std(axis=0)
   
    _parameter_std_dict = OrderedDict()
    for pn in _analyzer.parameter_names:
        _parameter_std_dict[pn] =_param_std_df.to_dict()[pn]
    
    return _parameter_std_dict

# testing methods to be developed

def test__get_parameter_variance__w_filenames():
    pass

def test__get_parameter_variance__w_objects():
    pass

def test__get_parameter_variance__w_filenames():
    pass

def test__get_parameter_variance__w_objects():
    pass

if __name__ == "__main__":
    _pypospack_root = pypospack.utils.get_pypospack_root_directory()
    _data_in_directory = os.path.join(_pypospack_root,'examples',
            'Ni__eam__born_exp_fs__sensitivityanalysis',
            'data__from_pareto_optimization')
    _pyposmat_data_fn = os.path.join(_data_in_directory,'pyposmat.kde.6.out')
    _pyposmat_config_fn = os.path.join(_data_in_directory,'pyposmat.config.in')

    _best_params = get_best_parameterization(
            config_fn=_pyposmat_config_fn,
            data_fn=_pyposmat_data_fn)
    _std_params = get_parameter_variance(
            config_fn=_pyposmat_config_fn,
            data_fn=_pyposmat_data_fn,
            metric_name='d_metric',
            n=100)

    # make new configuration file from the old configuration file

    from pypospack.pyposmat.data import PyposmatConfigurationFile
    _o_config = PyposmatConfigurationFile()
    _o_config.read(filename=_pyposmat_config_fn)

    pn1 = 'p_NiNi_phi0'
    pn2 = 'e_Ni_F0'

    # change to normal distribution for parameters of interest
    for pn in [pn1,pn2]:
        _o_config.configuration['sampling_dist'][pn] = ['normal',{'mu':_best_params[pn],'sigma':_std_params[pn]}]

    # change the rest of the free parameters to static
    for pn in _o_config.free_parameter_names:
        if pn != pn1 and pn != pn2:
            _o_config.configuration['sampling_dist'][pn] = ['equals',_best_params[pn]]

    _pyposmat_config_fn = os.path.join('data','pyposmat.config.in')

    # change the first iteration type to parametric
    _o_config.configuration['sampling_type'][0]['type'] = 'parametric'
    print(_o_config.configuration['sampling_type'])
    print("new pyposmat configuration file in:{}".format(_pyposmat_config_fn))
    _o_config.write(filename=_pyposmat_config_fn) 

