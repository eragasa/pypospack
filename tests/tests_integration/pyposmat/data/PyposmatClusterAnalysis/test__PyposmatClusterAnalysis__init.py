import pytest
import os

from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDatafile

from cluster_analysis import PyposmatClusterAnalysis

pyposmat_src_dir = os.path.join(
    '..','..','..','..',
    'data_test','Ni__eam__born_exp_fs_00','data__Ni__eam__born_exp_fs_03'
)
pyposmat_configuration_fn = os.path.join(
    pyposmat_src_dir,'pyposmat.config.in'
)
pyposmat_data_fn = os.path.join(
    pyposmat_src_dir,'pyposmat.kde.5.out'
)

def test__init():
    o = PyposmatClusterAnalysis()
    assert type(o) is PyposmatClusterAnalysis
    assert o.configuration_fn is None
    assert o.configuration is None
    assert o.data_fn is None
    assert o.data is None

def test__read_configuration_file():
    o = PyposmatClusterAnalysis()
    o.read_configuration(filename=pyposmat_configuration_fn)
    assert o.configuration_fn is pyposmat_configuration_fn
    assert o.configuration is PyposmatConfigurationFile

def test__read_data_file():
    o = PyposmatClusterAnalysis()
    o.read_datafile(filename=pyposmat_data_fn)
    assert o.data_fn is pyposmat_data_fn
    assert o.data is PyposmatDatafile
