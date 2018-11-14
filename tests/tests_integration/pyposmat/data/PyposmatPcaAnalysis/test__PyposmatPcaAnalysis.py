import pytest
import os
import pandas as pd
import pypospack.utils
from pypospack.pyposmat.data.pca_analysis import PyposmatPcaAnalysis

pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
pyposmat_config_fn = os.path.join(pypospack_root_dir,'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.config.in')
pyposmat_data_fn = os.path.join(pypospack_root_dir,'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.kde.4.out')

def test__init():
    o_pca = PyposmatPcaAnalysis()

def test__read_configuration():
    o_pca = PyposmatPcaAnalysis()
    o_pca.read_configuration(filename=pyposmat_config_fn)

def test__read_data():
    o_pca = PyposmatPcaAnalysis()
    o_pca.read_data(filename=pyposmat_data_fn)

    assert type(o_pca.data.df) is pd.DataFrame

def test__select_data():
    o_pca = PyposmatPcaAnalysis()
    o_pca.read_configuration(filename=pyposmat_config_fn)
    o_pca.read_data(filename=pyposmat_data_fn)
    o_pca.select_data(types=['qoi'])
    o_pca.select_data(types=['param'])
    o_pca.select_data(types=['qoi','param'])

    assert type(o_pca.df) is pd.DataFrame

def test__normalize_data():
    o_pca = PyposmatPcaAnalysis()
    o_pca.read_configuration(filename=pyposmat_config_fn)
    o_pca.read_data(filename=pyposmat_data_fn)
    o_pca.select_data(types=['qoi','param'])

