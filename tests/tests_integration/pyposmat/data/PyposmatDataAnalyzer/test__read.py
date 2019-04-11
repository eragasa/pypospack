import pytest
import os
from collections import OrderedDict

import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataAnalyzer

def get_testing_set():
    pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
    testing_set = OrderedDict()
    testing_set['data_directory'] = os.path.join(
            pypospack_root_dir,'examples','Si__sw',
            'dev__pareto_optimization','data')
    testing_set['config_fn'] = os.path.join(
            testing_set['data_directory'],'pyposmat.config.in')
    testing_set['badparameters_fn'] = os.path.join(
            testing_set['data_directory'],'pyposmat.badparameters.out')
    return testing_set

def dev__read():
    testing_set = get_testing_set()

    o = PyposmatDataAnalyzer(config_fn=
