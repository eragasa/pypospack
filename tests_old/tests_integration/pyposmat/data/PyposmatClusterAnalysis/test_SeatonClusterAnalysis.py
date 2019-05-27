import os
import pandas as pd
import pypospack.utils
from pypospack.pyposmat.data.cluster_analysis import SeatonClusterAnalysis

pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
pyposmat_config_fn = os.path.join(pypospack_root_dir, 'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.config.in')
pyposmat_data_fn = os.path.join(pypospack_root_dir, 'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.kde.4.out')


def test__init():
    o_cluster = SeatonClusterAnalysis()


def test__calculate_clusters():
    o_cluster = SeatonClusterAnalysis()
    o_cluster.read_configuration(filename=pyposmat_config_fn)
    o_cluster.read_data(filename=pyposmat_data_fn)
    o_cluster.calculate_clusters()
    assert type(o_cluster.df) is pd.DataFrame
