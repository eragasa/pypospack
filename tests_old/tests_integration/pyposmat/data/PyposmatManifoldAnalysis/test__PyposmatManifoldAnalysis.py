import os
import pandas as pd
import pypospack.utils
from pypospack.pyposmat.data.manifold_analysis import PyposmatManifoldAnalysis

pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
pyposmat_config_fn = os.path.join(pypospack_root_dir, 'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.config.in')
pyposmat_data_fn = os.path.join(pypospack_root_dir, 'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.kde.4.out')


def test__init():
    o_man = PyposmatManifoldAnalysis()


def test_new_transform_pca():
    o_man = PyposmatManifoldAnalysis()
    o_man.read_configuration(filename=pyposmat_config_fn)
    o_man.read_data(filename=pyposmat_data_fn)
    o_man.calculate_manifold()
    assert type(o_man.df) is pd.DataFrame
