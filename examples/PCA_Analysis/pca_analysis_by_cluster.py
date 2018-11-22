import os
import pandas as pd
import pypospack.utils
from pypospack.pyposmat.data.pipeline import PyposmatPipeline

pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
pyposmat_config_fn = os.path.join(pypospack_root_dir, 'examples/PCA_Analysis/configure_param_clustering.in')
pyposmat_data_fn = os.path.join(pypospack_root_dir, 'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.kde.4.out')


if __name__ == "__main__":
    pipeline = PyposmatPipeline(configuration_fn=pyposmat_config_fn,
                                data_fn=pyposmat_data_fn)
    pipeline.run()
    print(pipeline.df)
