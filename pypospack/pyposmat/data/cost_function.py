import pandas as pd

from collections import PyposmatDataFile

class CostFunction(object):

    def __init__(self,configuration_fn=None,datafile_fn=None):

        if configuration_fn is not None:
            self.read_configuration_file(filename=configuration_fn)
        else:
            self.configuration_fn = None
            self.configuration = None

        if datafile_fn is not None:
            self.datafile_fn = datafile_fn
            self.datafile = PyposmatDataFile()
            self.datafie.read(filename=self.datafile_fn)
        else:
            self.datafile_fn = None
            self.datafile = None

        self.weights = None
        self.qoi_targets = None
    
    def read_configuration_file(self,filename):
        self.configuration_fn = filename
        self.configuration = PyposmatConfigurationFile()
        self.configuration.read(filename=self.configuration_fn)

    def set_weights(self,weights):

        self.weights = weights

    def set_qoi_targets(self,qoi_targets):

        self.qoi_targets(qoi_targets)

    def evaluate_with_dataframe(self,df):

        if not isintance(df,pd.DataFrame):
            raise ValueError('df must be a pandas.DataFrame')

    def evaluate_with_datafile(self,datafile):

        if type(datafile) is str:
            _datafile = PyposmatDataFile()
            _datafile.read(datafile)
        elif isinstance(datafile,PyposmatDataFile):
            pass
