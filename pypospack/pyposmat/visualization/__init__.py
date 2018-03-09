
import copy,os
from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib import cm

from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

class PyposmatDatafileVisualization(object):
    def __init__(self):
        self.datafile = None
        self.configuation = None
    @property
    def parameter_names(self):
        return self._parameter_names

    @property
    def qoi_names(self):
        return self._qoi_names

    @property
    def error_names(self):
        return self._error_names

    @property
    def df(self):
        return self._df

    @property
    def parameter_df(self):
        return self._parameter_df

    @property
    def qoi_df(self):
        return self._qoi_df

    @property
    def error_df(self):
        error_df(self)
        return self._error_df

    def read_configuration(self,filename):
        self.configuration = PyposmatConfigurationFile()
        self.configuration.read(filename=filename)
    
    def read_datafile(self,filename):
        self.datafile = PyposmatDataFile()
        self.datafile.read(filename=filename)

        self._parameter_names = self.datafile.parameter_names
        self._qoi_names = self.datafile.qoi_names
        self._error_names = self.datafile.error_names

        self._df = copy.deepcopy(self.datafile.df)
        self.create_absolute_errors()

    def create_absolute_errors(self):
        for q in self.qoi_names:
            aen = "{}.abserr".format(q)
            en = "{}.err".format(q)
            self._df[aen] = self._df[en].abs()

