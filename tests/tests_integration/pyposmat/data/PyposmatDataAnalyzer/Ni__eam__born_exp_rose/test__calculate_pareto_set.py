import os
import pytest

import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataAnalyzer

config_fn = os.path.join(
        pypospack.utils.get_pypospack_root_directory(),
        'examples','Ni__eam__born_exp_rose','01_preconditioning_3.5NN','data',
        'pyposmat.config.in')

data_fn = os.path.join(
        pypospack.utils.get_pypospack_root_directory(),
        'examples','Ni__eam__born_exp_rose','01_preconditioning_3.5NN','data',
        'pyposmat.results.0.out')

def test__calculate_pareto_set__df_equals_DataFrame():
    o = PyposmatDataAnalyzer(fn_config=config_fn,fn_data=data_fn)
    o.calculate_pareto_set(df=o.data.df)

def test__calculate_pareto_set__df_equals_None():
    o = PyposmatDataAnalyzer(fn_config=config_fn,fn_data=data_fn)
    o.calculate_pareto_set()

