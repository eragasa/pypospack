import os
import pandas as pd
import pypospack.utils
from pypospack.pyposmat.data.pca_analysis import PyposmatPcaAnalysis

pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
pyposmat_config_fn = os.path.join(pypospack_root_dir, 'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.config.in')
pyposmat_data_fn = os.path.join(pypospack_root_dir, 'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.kde.4.out')


def test__init():
    o_pca = PyposmatPcaAnalysis()


def test_new_transform_pca():
    o_pca = PyposmatPcaAnalysis()
    o_pca.read_configuration(filename=pyposmat_config_fn)
    o_pca.read_data(filename=pyposmat_data_fn)
    o_pca.transform_pca()
    assert type(o_pca.df) is pd.DataFrame
