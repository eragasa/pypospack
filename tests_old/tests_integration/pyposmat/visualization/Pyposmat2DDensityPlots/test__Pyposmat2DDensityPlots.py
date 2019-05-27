import pytest
import matplotlib.pyplot as plt
import os,shutil
import pypospack.utils
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.visualization import Pyposmat2DDensityPlot

config_fn = os.path.join(
        pypospack.utils.get_pypospack_root_directory(),
        'data','Si__sw__data',
        'pareto_optimization_unconstrained',
        'pyposmat.config.in')
data_fn = os.path.join(
        pypospack.utils.get_pypospack_root_directory(),
        'data','Si__sw__data',
        'pareto_optimization_unconstrained',
        'pyposmat.kde.20.out')

assert os.path.isfile(config_fn)
assert os.path.isfile(data_fn)

def test____init____no_args():
    o = Pyposmat2DDensityPlot()

    assert isinstance(o,Pyposmat2DDensityPlot)
    assert o.configuration is None
    assert o.data is None
    assert o.fig is None
    assert o.ax is None

def test____init____config_as_obj():
    o = Pyposmat2DDensityPlot(config=PyposmatConfigurationFile())

    assert isinstance(o,Pyposmat2DDensityPlot)
    assert isinstance(o.configuration,PyposmatConfigurationFile)
    assert o.data is None
    assert o.fig is None
    assert o.ax is None

def test____init____config_as_str():
    o = Pyposmat2DDensityPlot(config=config_fn)

    assert isinstance(o,Pyposmat2DDensityPlot)
    assert isinstance(o.configuration,PyposmatConfigurationFile)
    assert o.data is None
    assert o.fig is None
    assert o.ax is None

def test____init____data_as_obj():
    o = Pyposmat2DDensityPlot(data=PyposmatDataFile())

    assert isinstance(o,Pyposmat2DDensityPlot)
    assert o.configuration is None
    assert isinstance(o.data,PyposmatDataFile)
    assert o.fig is None
    assert o.ax is None

def test____init____data_as_str():
    o = Pyposmat2DDensityPlot(data=data_fn)

    assert isinstance(o,Pyposmat2DDensityPlot)
    assert o.configuration is None
    assert isinstance(o.data,PyposmatDataFile)
    assert o.fig is None
    assert o.ax is None

def test____init____all_args_as_obj():
    o = Pyposmat2DDensityPlot(config=PyposmatConfigurationFile(),
                              data=PyposmatDataFile())

    assert isinstance(o,Pyposmat2DDensityPlot)
    assert isinstance(o.configuration,PyposmatConfigurationFile)
    assert isinstance(o.data,PyposmatDataFile)
    assert o.fig is None
    assert o.ax is None

def test____init____all_args_as_str():
    o = Pyposmat2DDensityPlot(config=config_fn,
                              data=data_fn)

    assert isinstance(o,Pyposmat2DDensityPlot)
    assert isinstance(o.configuration,PyposmatConfigurationFile)
    assert isinstance(o.data,PyposmatDataFile)
    assert o.fig is None
    assert o.ax is None

def test__create_subplots__no_args():
    o = Pyposmat2DDensityPlot(config=config_fn,data=data_fn)
    o.create_subplots()

    assert isinstance(o.fig,plt.Figure)
    assert isinstance(o.ax,plt.Axes)

def dev__determine_limits():
    o = Pyposmat2DDensityPlot(config=config_fn,data=data_fn)
    o.create_subplots()
    for qoi_name in o.configuration.qoi_names:
        print(qoi_name)
        print(o.determine_limits(qoi_name))

def test__determine_limits():
    o = Pyposmat2DDensityPlot(config=config_fn,data=data_fn)
    for qoi_name in o.configuration.qoi_names:
        lim_min,lim_max = o.determine_limits(qoi_name)
        assert isinstance(lim_min,float)
        assert isinstance(lim_max,float)
        assert lim_min < lim_max

def dev__determine_x_limits():
    o = Pyposmat2DDensityPlot(config=config_fn,data=data_fn)
    o.x_name = o.configuration.qoi_names[0]

    o.determine_x_limits()
    print(o.x_limits)

def test__determine_x_limits():
    o = Pyposmat2DDensityPlot(config=config_fn,data=data_fn)
    o.x_name = o.configuration.qoi_names[0]

    o.determine_x_limits()

    assert isinstance(o.x_limits[0],float)
    assert isinstance(o.x_limits[1],float)

def dev__determine_y_limits():
    o = Pyposmat2DDensityPlot(config=config_fn,data=data_fn)
    o.y_name = o.configuration.qoi_names[1]

    o.determine_y_limits()
    print(o.y_limits)

def test__determine_y_limits():
    o = Pyposmat2DDensityPlot(config=config_fn,data=data_fn)
    o.y_name = o.configuration.qoi_names[1]

    o.determine_y_limits()

    assert isinstance(o.y_limits[0],float)
    assert isinstance(o.y_limits[1],float)

def test__plot():
    o = Pyposmat2DDensityPlot(config=config_fn,data=data_fn)
    x_name = o.configuration.qoi_names[0]
    y_name = o.configuration.qoi_names[1]
    o.plot(x_name=x_name,y_name=y_name)
    plt.show()

if __name__ == "__main__":
    dev__determine_x_limits()
    dev__determine_y_limits()
