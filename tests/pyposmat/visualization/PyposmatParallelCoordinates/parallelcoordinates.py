import copy, os
import numpy as np
import pandas as pd
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

class PyposmatParallelCoordinates(object):
    
    def __init__(self):
        self._configuration = PyposmatConfigurationFile()
        self._data = PyposmatDataFile()

    def set_configuration(self, configuration):
        if isinstance(configuration, str):
            self.set_configuration_by_path(path=configuration)
        elif isinstance(configuration, PyposmatConfigurationFile):
            self.set_configuration_by_object(config_obj=configuration)
        else:
            raise TypeError

    def set_data(self, data):
        if isinstance(data, str):
            self.data = PyposmatDataFile()
            self.data.read(data)
        elif isinstance(data, PyposmatDataFile):
            self.data = data
        else:
            raise TypeError
    
    def set_configuration_by_path(self, path):
        assert isinstance(path, str):
        self.configuration = PyposmatConfigurationFile()
        self.configuration.read(path)

    def set_configuration_by_obj(self, config_obj):
        assert isinstance(configuration, PyposmatConfigurationFile)
        self.configuration = config_obj

    def set_data_by_path(self, path):
        assert isinstance(path, str)
        self.data  = PyposmatDataFile()
        self.data.read(path)

    def set_data_by_obj(self, data_obj):
        assert isinstance(data_obj, PyposmatDataFile)
        self.data = data

    @property
    def configuration(self):
        return self._configuration

    @configuration.setter
    def configuration(self, configuration):
        assert isinstance(configuration, PyposmatConfigurationFile)
        self._configuration = configuration

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, data):
        assert isinstance(configuration, PyposmatDataFile)
        self._data = data

    def plot(self, path):
        assert isinstance(path, str)

