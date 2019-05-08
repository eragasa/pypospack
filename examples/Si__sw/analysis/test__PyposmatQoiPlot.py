import pytest

from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile

from abstract_plot import PyposmatAbstractPlot
from qoi_plot import PyposmatQoiPlot



def test__init__no_args():
    o_plot = PyposmatQoiPlot()
    
    assert isinstance(o_plot,PyposmatAbstractPlot)

def test__init__config_as_obj():
    o_plot = PyposmatQoiPlot(config=PyposmatConfigurationFile())

    assert isinstance(o_plot.configuration,PyposmatConfigurationFile)
