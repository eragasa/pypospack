import pytest,os
import numpy as np
import pandas as pd
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile

pyposmat_root = [v for v in os.environ['PYTHONPATH'].strip().split(':') if v.endswith('pypospack')][0]

class PyposmatPostProcessorTestHarness(object):
    
    def __init__(self,configuration_fn,datafile_fn):
        self.configuration_fn = configuration_fn
        self.datafile_fn = datafile_fn

        if configuration_fn is not None:
            self.configuration = PyposmatConfigurationFile()
            self.configuration.read(configuration_fn)
        if datafile_fn is not None:
            self.datafile = PyposmatDataFile()
            self.datafile.read(filename=datafile_fn)

    def get_parameter_names(self):
        return self.configuration.parameter_names

data_directory = os.path.join(pyposmat_root,'data','MgO_pareto_data')
config_fn = os.path.join(data_directory,'pyposmat.config.in')
data_fn = os.path.join(data_directory,"culled_004.out")

from post_processor import PyposmatPostProcessor

def test__harness():
    assert os.path.isfile(config_fn)
    assert os.path.isfile(data_fn)

def test____init__():
    o = PyposmatPostProcessor()

    assert o._configuration is None
    assert o._datafile is None
    assert o.configuration is None
    assert o.datafile is None

def test____init__with_config_fn():
    o = PyposmatPostProcessor(configuration_fn=config_fn)

    assert type(o.configuration) is PyposmatConfigurationFile
    assert o.datafile is None

def test____init__with_data_fn():
    o = PyposmatPostProcessor(datafile_fn=data_fn)

    assert o.configuration is None
    assert type(o.datafile) is PyposmatDataFile
    assert o.datafile_fn == data_fn

    assert isinstance(o.df,pd.DataFrame)
def test____init__with_config_fn_and_data_fn():
    o = PyposmatPostProcessor(
            configuration_fn=config_fn,
            datafile_fn=data_fn)

    assert type(o.configuration) is PyposmatConfigurationFile
    assert type(o.datafile) is PyposmatDataFile
    assert o.configuration_fn == config_fn
    assert o.datafile_fn == data_fn


def test__read_configuration():
    o = PyposmatPostProcessor()
    o.read_configuration(filename=config_fn)
    
    assert type(o.configuration) is PyposmatConfigurationFile
    assert o.configuration_fn == config_fn
    
    assert o.datafile is None
    assert o.datafile_fn == None

    assert type(o.qoi_names) is list
    assert type(o.error_names) is list
    assert type(o.qoi_validation_names) is list
    assert type(o.error_validation_names) is list
    assert type(o.qoi_testing_names) is list
    assert type(o.error_testing_names) is list
    assert type(o.qoi_fitting_names) is list
    assert type(o.error_testing_names) is list

def test__read_datafile():
    o = PyposmatPostProcessor()
    o.read_datafile(filename=data_fn)

    assert o.configuration is None
    assert o.configuration_fn is None

    assert type(o.datafile) is PyposmatDataFile
    assert o.datafile_fn == data_fn
if __name__ == "__main__":
    o = PyposmatPostProcessor()
