import pytest

from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile
from abstract_plot import PyposmatAbstractPlot

def test__init__no_args():
    o_plot = PyposmatAbstractPlot(config=None,data=None)
    assert o_plot.configuration is None
    assert o_plot.data is None
    assert o_plot.fig is None
    assert o_plot.ax is None

def test__init__none_args():
    o_plot = PyposmatAbstractPlot()
    assert o_plot.configuration is None
    assert o_plot.data is None
    assert o_plot.fig is None
    assert o_plot.ax is None

def test__init__config_is_obj():
    o_plot = PyposmatAbstractPlot(config=PyposmatConfigurationFile())
    assert isinstance(o_plot.configuration,PyposmatConfigurationFile)
    assert o_plot.data is None
    assert o_plot.fig is None
    assert o_plot.ax is None

def test__init__data_is_obj():
    o_plot = PyposmatAbstractPlot(data=PyposmatDataFile())
    assert o_plot.configuration is None
    assert isinstance(o_plot.data,PyposmatDataFile)
    assert o_plot.fig is None
    assert o_plot.ax is None

def dev__init__no_args_1():
    o_plot = PyposmatAbstractPlot(config=None,data=None)
    assert o_plot.configuration is None
    assert o_plot.data is None
    assert o_plot.fig is None
    assert o_plot.ax is None



