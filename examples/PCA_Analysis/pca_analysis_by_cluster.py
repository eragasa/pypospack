import os
import pandas as pd
import pypospack.utils
from pypospack.pyposmat.data.cluster_analysis import SeatonClusterAnalysis
from pypospack.pyposmat.data.preprocess import PyposmatPreprocessor
from pypospack.pyposmat.data.pca_analysis import PyposmatPcaAnalysis

pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
pyposmat_config_fn = os.path.join(pypospack_root_dir, 'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.config.in')
pyposmat_data_fn = os.path.join(pypospack_root_dir, 'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.kde.4.out')


if __name__ == "__main__":
    o_pre = PyposmatPreprocessor()
    o_pre.read_configuration(filename=pyposmat_config_fn)
    o_pre.read_data(filename=pyposmat_data_fn)

    o_pca = PyposmatPcaAnalysis()
    o_pca.read_configuration(filename=pyposmat_config_fn)
    o_pca.read_data(filename=pyposmat_data_fn)

    o_cluster = SeatonClusterAnalysis()
    o_cluster.read_configuration(filename=pyposmat_config_fn)
    o_cluster.read_data(filename=pyposmat_data_fn)

    o_pre.select_data()
    norm_df = o_pre.normalize_data()
    pca_df = o_pca.transform_pca(norm_df)
    tsne_df = o_cluster.calculate_manifold(pca_df)
    cluster_df = o_cluster.calculate_clusters(tsne_df)

    df = pd.concat([pca_df, cluster_df])
    print(df.shape)

