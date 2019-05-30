import pytest

from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile

from pypospack.pyopsmat.visualization import PyposmatAbstractPlot
from pypospack.pyposmat.visualization import PyposmatQoiPlot

def test__init__no_args():
    o_plot = PyposmatQoiPlot()
    
    assert isinstance(o_plot,PyposmatAbstractPlot)

def test__init__config_as_obj():
    o_plot = PyposmatQoiPlot(config=PyposmatConfigurationFile())

    assert isinstance(o_plot.configuration,PyposmatConfigurationFile)
