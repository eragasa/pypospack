import os
import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.visualization import PyposmatAbstractPlot
from parallel_coordinates_plot import PyposmatParallelCoordinatesPlot

import matplotlib.pyplot as plt
def test____init____no_args():
    o = PyposmatParallelCoordinatesPlot()
    assert isinstance(o,PyposmatAbstractPlot)
    assert o.configuration is None
    assert o.data is None
    assert o.excluded_names == []

    assert o.fig is None
    assert o.ax is None

def test__create_subplots():
    o = PyposmatParallelCoordinatesPlot()
    o.create_subplots()

    assert isinstance(o.fig,plt.Figure)
    assert isinstance(o.ax,plt.Axes)

def test__initalize_configuration__config_None():
    o = PyposmatParallelCoordinatesPlot()
    o.initialize_configuration(config=None)
    assert o.configuration is None

def test__initialize_configuration__config_obj():
    o = PyposmatParallelCoordinatesPlot()
    o.initialize_configuration(config=PyposmatConfigurationFile())
    assert isinstance(o.configuration,PyposmatConfigurationFile)

def test__initalize_configuration__config_str():
    config_fn = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'data','Si__sw__data',
            'pareto_optimization_unconstrained')

    o = PyposmatParallelCoordinatesPlot()
    o.initial_configuration(config=config_fn)
    assert isinstance(o.configuration,PyposmatConfigurationFile)

def dev__create_subplots():
    o = PyposmatParallelCoordinatesPlot()
    o.create_subplots()

    print(type(o.fig))
    print(type(o.ax))

if __name__ == "__main__":
    dev__create_subplots()
