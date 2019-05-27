import pytest
import os
from collections import OrderedDict
import pandas as pd
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

def dev__read():
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

    o.read(filename=testing_set['badparameters_in_fn'])

    print(type(o.df))
    print(isinstance(o.df,pd.DataFrame))
    print(o.df)

def test__read():
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

    o.read(filename=testing_set['badparameters_in_fn'])

    assert isinstance(o.df,pd.DataFrame)
if __name__ == "__main__":
    dev__read()
