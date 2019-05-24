from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
class PyposmatDataFileVisualization(object):
    def __init__(self):
        self._datafile = None
        self._configuation = None
    
    @property
    def configuration(self):return self._configuration
    @configuration.setter
    def configuration(self,config):
        self._configuration = config
    @property
    def datafile(self):return self._datafile
    @datafile.setter
    def datafile(self,datafile):
        self._datafile=datafile
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

    def read_datafiles(self,
            filenames):

        if type(filenames) is not list:
            raise ValueError("filenames must be a list")
        self.n_iterations = len(filenames)
    
        self.datafiles = OrderedDict()
        self._df = OrderedDict()
        for k in filenames:
            self.datafiles[k] = PyposmatDataFile()
            self.datafiles[k].read(k)
            self._parameter_names = self.datafiles[k].parameter_names
            self._qoi_names = self.datafiles[k].qoi_names
            self._error_names = self.datafiles[k].error_names
            self._df[k] = copy.deepcopy(self.datafiles[k].df)

            for q in self._qoi_names:
                aen = "{}.abserr".format(q)
                en = "{}.err".format(q)
                self._df[k][aen] = self._df[k][en].abs()
    
    def create_absolute_errors(self):
        for q in self.qoi_names:
            aen = "{}.abserr".format(q)
            en = "{}.err".format(q)
            self._df[aen] = self._df[en].abs()

