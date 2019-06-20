# author: Eugene J. Ragasa <eragasa@ufl.edu>
# License: MIT
# Copyright: 2019

from abc import ABC, abstractmethod

class AbstractAnalysis(ABC):
    def __init__(self,
                 pyposmat_configuration,
                 pyposmat_data,
                 names):
        super().__init__()
        self.configuration = None
        self.data = None
        self.names = None
        self.initialize_configuration(configuration=pyposmat_configuration)
        self.initialize_data(data=pyposmat_data)
        self.initialize_names(names=names)
    
    def initialize_configuration(self,configuration):
        if isinstance(configuration,PyposmatConfigurationFile):
            self.configuration = configuration
        elif isinstance(configuration, str):
            self.configuration = PyposmatConfigurationFile()
            self.configuration.read(filename=configuration)
        else:
            msg = ("configuration must either be a path to a configuration "
                   "file or an instance of PyposmatConfigurationFile")
            raise TypeError(msg)

    def initialize_data(self, data):
        if isinstance(data, PyposmatDataFile):
            self.data = data
        elif isinstance(data, str):
            self.data = PyposmatDataFile()
            self.data.read(filename=data)
        else:
            msg = ("data must either be a path to a data file or an instance "
                    "of PyposmatDataFile")
            raise TypeError(msg)

    def initialize_names(self,names):
        if isinstance(names,list):
            self.names = names
        elif isinstance(names,str):
            if names == 'all':
                self.names = self.configuration.free_parameter_names\
                        + self.configuration.qoi_names
            elif names == 'free_parameters':
                self.names = self.configuration.free_parameter_names
            elif names == 'qois':
                self.names = self.configuration.qoi_names
            else:
                msg = ("valid names arguments are all,free_parameters,qois or "
                       "a list of names")
                raise ValueError(msg)
        else:
            msg = ("valid names arguments are all,free_parameters,qois or "
                   "a list of names")
            raise TypeError(msg)

