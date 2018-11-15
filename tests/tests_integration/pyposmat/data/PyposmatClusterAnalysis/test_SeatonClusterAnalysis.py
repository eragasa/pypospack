import os
import pandas as pd
import pypospack.utils
from pypospack.pyposmat.data.cluster_analysis import SeatonClusterAnalysis
from pypospack.pyposmat.data.preprocess import PyposmatPreprocessor
from pypospack.pyposmat.data.pca_analysis import PyposmatPcaAnalysis

pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
pyposmat_config_fn = os.path.join(pypospack_root_dir, 'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.config.in')
pyposmat_data_fn = os.path.join(pypospack_root_dir, 'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.kde.4.out')


def test__init():
    o_cluster = SeatonClusterAnalysis()


def test__read_configuration():
    o_cluster = SeatonClusterAnalysis()
    o_cluster.read_configuration(filename=pyposmat_config_fn)


def test__read_data():
    o_cluster = SeatonClusterAnalysis()
    o_cluster.read_data(filename=pyposmat_data_fn)
    assert type(o_cluster.data.df) is pd.DataFrame


def test__calculate_manifold():
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

    o_cluster = SeatonClusterAnalysis()
    o_cluster.read_configuration(filename=pyposmat_config_fn)
    o_cluster.read_data(filename=pyposmat_data_fn)
    manifold_df = o_cluster.calculate_manifold(pca_df)  # calculate tsne manifold
    assert type(manifold_df) is pd.DataFrame

    cluster_df = o_cluster.calculate_clusters(manifold_df)  # calculate cluster labels in tsne space
    assert type(cluster_df) is pd.DataFrame

    cluster_0_df = o_cluster.select_cluster(0, cluster_df)  # select cluster id 0 from larger df
    assert type(cluster_0_df) is pd.DataFrame
