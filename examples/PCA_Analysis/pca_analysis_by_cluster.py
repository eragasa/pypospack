import os
import pandas as pd
import pypospack.utils
from pypospack.pyposmat.data.pipeline import PyposmatPipeline

pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
pyposmat_param_config_fn = os.path.join(pypospack_root_dir, 'examples/PCA_Analysis/configure_param_clustering.in')
pyposmat_qoi_config_fn = os.path.join(pypospack_root_dir, 'examples/PCA_Analysis/configure_qoi_clustering.in')

pyposmat_data_fn = os.path.join(pypospack_root_dir, 'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.kde.4.out')


if __name__ == "__main__":
    # cluster the parameter space with kmeans
    param_pipeline = PyposmatPipeline(configuration_fn=pyposmat_param_config_fn,
                                      data_fn=pyposmat_data_fn)
    param_pipeline.read_configuration(pyposmat_param_config_fn)
    param_pipeline.read_data(pyposmat_data_fn)
    param_pipeline.run()

    # pca transform qoi space
    qoi_pipeline = PyposmatPipeline(configuration_fn=pyposmat_qoi_config_fn,
                                    data_fn=pyposmat_data_fn)
    qoi_pipeline.read_configuration(pyposmat_qoi_config_fn)
    qoi_pipeline.read_data(pyposmat_data_fn)
    qoi_pipeline.run()

    qoi_pipeline.df['cluster_id'] = param_pipeline.df['cluster_id']  # transfer cluster labels
