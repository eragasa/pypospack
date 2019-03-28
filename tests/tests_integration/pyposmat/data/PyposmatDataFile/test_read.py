import pytest
import os
from collections import OrderedDict

import pypospack.utils
from pypospack.pyposmat import PyposmatDataFile

def get_testing_set():
    testing_set = OrderedDict()
    testing_set['results_data_fn'] = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'data','Si__sw__data','pareto_optimization_1','data','pyposmat.results.0.out'
    )

    assert os.path.isfile(testing_set['results_data_fn'])

    return testing_set

def dev__read():
    
    testing_set = get_testing_set()

    o = PyposmatDataFile()
    o.read

