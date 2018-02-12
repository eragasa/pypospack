import os
from pypospack.pyposmat import PyposmatConfigurationFile
from pypospack.pyposmat import PyposmatDataFile
import pareto


class PypospackDataAnalyzer(object):
    def __init__(self):
        self._pyposmat_configuration = None

    @property
    def pyposmat_configuration(self):
        return self._pypospack_configuration

    @pyposmat_configuration.setter
    def pypospack_configuration(self,configuration):
        assert type(configuration) is PypospackConfigurationFile


    def read_pypospack_configuration_file(self,filename):
        self._pyposmat_configuation = PyposmatConfigurationFile()

if __name__ == "__main__":
    data_directory = 'data'
    pyposmat_data_filename = 'pypospack.results.out'
   
    datafile = PyposmatDataFile(filename=os.path.join(
        data_directory,pyposmat_data_filename))

    datafile.read()

    print(datafile.parameter_names)
    print(datafile.qoi_names)
    print(datafile.error_names)
    print(datafile.df)
 
