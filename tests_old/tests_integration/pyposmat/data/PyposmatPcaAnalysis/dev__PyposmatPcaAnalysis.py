import pytest
import pypospack.utils


pypospack_root_dir = get_pypospack_root_directory()
pyposmat_config_file = os.path.join(pypospack_root_dir,'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.config.in')
pyposmat_data_file = os.path.join(pypospack_root_dir,'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.results.4.out')

def test__init():
    o_pca = PyposmatPcaAnalysis()


