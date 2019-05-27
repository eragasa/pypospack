import os
import pandas as pd
import pypospack.utils
from pypospack.pyposmat.data.pipeline import BasePipeSegment
from pypospack.pyposmat.data.pipeline import PyposmatPipeline

pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
pyposmat_config_fn = os.path.join(pypospack_root_dir, 'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.config.in')
pyposmat_data_fn = os.path.join(pypospack_root_dir, 'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.kde.4.out')


def test__init():
    o_bps = BasePipeSegment()


def test__read_configuration():
    o_bps = BasePipeSegment()
    o_bps.read_configuration(filename=pyposmat_config_fn)


def test__read_data():
    o_bps = BasePipeSegment()
    o_bps.read_data(filename=pyposmat_data_fn)
    assert type(o_bps.data.df) is pd.DataFrame


def test__select_data():
    o_bps = BasePipeSegment()
    o_bps.read_configuration(filename=pyposmat_config_fn)
    o_bps.read_data(filename=pyposmat_data_fn)
    o_bps.select_data()
    o_bps.select_data(types=['qoi'])
    o_bps.select_data(types=['err'])
    o_bps.select_data(types=['qoi', 'param'])
    assert type(o_bps.df) is pd.DataFrame
