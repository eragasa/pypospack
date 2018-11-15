import os

from collections import OrderedDict
import pandas as pd
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatClusterAnalysis


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
    assert type(o.configuration) is PyposmatConfigurationFile


def test__read_data():
    o = PyposmatClusterAnalysis()
    o.read_data(filename=pyposmat_data_fn)
    assert o.data_fn is pyposmat_data_fn
    assert type(o.data) is PyposmatDataFile
    assert type(o.parameter_names) is list
    assert type(o.qoi_names) is list
    assert type(o.error_names) is list
    assert type(o.df) is pd.DataFrame


def test__init_from_ordered_dict__minimum_config():
    d = OrderedDict([
        ('configuration_fn',pyposmat_configuration_fn),
        ('data_fn',pyposmat_data_fn)
    ])
    o = PyposmatClusterAnalysis.init_from_ordered_dict(d)

    # check configuration
    assert o.configuration_fn is pyposmat_configuration_fn
    assert type(o.configuration) is PyposmatConfigurationFile

    # check data
    assert o.data_fn is pyposmat_data_fn
    assert type(o.data) is PyposmatDataFile
    assert o.include_parameters is True
    assert o.include_qois is False
    assert o.include_errors is False


def test__init_from_orderedDict__variable_selection():
    d = OrderedDict([
        ('configuration_fn',pyposmat_configuration_fn),
        ('data_fn',pyposmat_data_fn),
        ('include_parameters',True),
        ('include_qois', False),
        ('include_errors',False)
    ])

    o = PyposmatClusterAnalysis.init_from_ordered_dict(d)

    # check configuration
    assert o.configuration_fn is pyposmat_configuration_fn
    assert type(o.configuration) is PyposmatConfigurationFile
    assert o.include_parameters is d['include_parameters']
    assert o.include_qois is d['include_qois']
    assert o.include_errors is d['include_errors']

    assert type(o.parameter_names) is list
    assert type(o.qoi_names) is list
    assert type(o.error_names) is list

    assert o.cluster_names == o.parameter_names
