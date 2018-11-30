import os
import pypospack.utils
from pypospack.pyposmat.data.pipeline import PyposmatPipeline

pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
configuration_dir = 'examples/Ni__eam__born_exp_fs_postprocessing/Reduced_TSNE_qoi_in_param/configuration/'

config_fn_0 = os.path.join(pypospack_root_dir,
                           configuration_dir,
                           'configure_qoi_clustering.in')
config_fn_1 = os.path.join(pypospack_root_dir,
                           configuration_dir,
                           'configure_param_tsne.in')
config_fn_2 = os.path.join(pypospack_root_dir,
                           configuration_dir,
                           'configure_param_plot.in')

pyposmat_data_fn = os.path.join(pypospack_root_dir,
                                'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.kde.4.out')


if __name__ == "__main__":
    # normalize, pca transform, and cluster qoi space
    pipeline_0 = PyposmatPipeline(configuration_fn=config_fn_0,
                                  data_fn=pyposmat_data_fn)
    pipeline_0.read_configuration(config_fn_0)
    pipeline_0.read_data(pyposmat_data_fn)
    pipeline_0.run()

    # normalize and pca transform param space
    pipeline_1 = PyposmatPipeline(configuration_fn=config_fn_1,
                                  data_fn=pyposmat_data_fn)
    pipeline_1.read_configuration(config_fn_1)
    pipeline_1.read_data(pyposmat_data_fn)
    pipeline_1.run()

    pipeline_1.df['cluster_id'] = pipeline_0.df['cluster_id']  # transfer cluster labels

    # plot parameter clusters in param space
    pipeline_2 = PyposmatPipeline(configuration_fn=config_fn_2,
                                  df=pipeline_1.df)
    pipeline_2.read_configuration(config_fn_2)
    pipeline_2.run()
