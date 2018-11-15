import os
import pandas as pd
import pypospack.utils
from pypospack.pyposmat.data.pca_analysis import PyposmatPcaAnalysis
from pypospack.pyposmat.data.preprocess import PyposmatPreprocessor

pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
pyposmat_config_fn = os.path.join(pypospack_root_dir, 'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.config.in')
pyposmat_data_fn = os.path.join(pypospack_root_dir, 'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.kde.4.out')


def test__init():
    o_pca = PyposmatPcaAnalysis()


def test__read_configuration():
    o_pca = PyposmatPcaAnalysis()
    o_pca.read_configuration(filename=pyposmat_config_fn)


def test__read_data():
    o_pca = PyposmatPcaAnalysis()
    o_pca.read_data(filename=pyposmat_data_fn)
    assert type(o_pca.data.df) is pd.DataFrame


def test_new_transform_pca():
    o_pre = PyposmatPreprocessor()
    o_pre.read_configuration(filename=pyposmat_config_fn)
    o_pre.read_data(filename=pyposmat_data_fn)
    o_pre.select_data()
    norm_df = o_pre.normalize_data()  # init with no ordered dict for default settings
    assert type(norm_df) is pd.DataFrame

    o_pca = PyposmatPcaAnalysis()
    o_pca.read_configuration(filename=pyposmat_config_fn)
    o_pca.read_data(filename=pyposmat_data_fn)
    pca_df = o_pca.transform_pca(norm_df)  # transform the normalized data
    assert type(pca_df) is pd.DataFrame
