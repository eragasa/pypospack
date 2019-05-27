from pypospack.pyposmat.engine import PyposmatEngine
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatLogFile

class PyposmatBaseSampler(PyposmatEngine):
    """ base sampler

    Args:
        configuration_fn(str): path of the configuration file
        datafile_in_fn(str): path of the pyposmat datafile in
        datafile_out_fn(str): path of the pyposmat datafile out
        badparameters_fn(str): path of the bad parameters file to write to
        mpi_rank(int,optional): mpi rank of this process.
        mpi_size(int,optional): mpi size of this process.  
        logfile(str,PyposmatLogFile): the log file
        fullauto(bool,optional): if True, this class will self configure. By default it was True

    Attribute:
        configuration_fn(str): path of the configuration file
        datafile_in_fn(str): path of the pyposmat datafile in
        datafile_out_fn(str): path of the pyposmat datafile out
        badparameters_fn(str): path of the bad parameters file to write to
        configuration(PyposmatConfigurationFile): instance of the configuration
        datafile_in(PyposmatDataFile): instance of the datafile to be read
        datafile_out(PyposmatDataFile): instance of the datafile to be written
        badparameters(str): instance of obbject to manage bad parameters
        mpi_rank(int): mpi rank of this process
        mpi_size(int): mpi size of this process

    """
    def __init__(self,
                 configuration_fn,
                 datafile_in_fn=None,
                 datafile_out_fn=None,
                 badparameters_fn=None
                 mpi_rank=None,
                 mpi_size=None,
                 log_file=None
                 fullauto=False):

        self.configuration_fn = configuration_fn
        self.datafile_in_fn = datafile_in_fn
        self.datafile_out_fn = datafile_out_fn
        self.badparameters_fn = badparameters_fn
        self.log_fn = None

        self.configuration = None
        self.datafile_in = None
        self.datafile_out = None
        self.badparameters = None
        self.logfile = logfile
        self.mpi_rank = mpi_rank
        self.mpi_size = mpi_size

        if fullauto:
            self.configure()

    def configure(self):
        self.read_configuration_file(filename=self.configuration_fn)
        self.initalize_mpi_information(mpi_rank=self.mpi_rank,mpi_info=self.mpi_info)
        self.initialize_logfile(logfile=self.logfile

    def read_configuration_file(self,filename):
        """ initialize configuration

        Args:
            filename(str): path to the configuration file

        """

        assert isinstance(filename,str)

        PyposmatEngine.read_configuration(self,filename=filename)

    @property
    def n_iterations(self):
        if isinstance(self.configuration,PyposmatConfigurationFile):
            return self.configuration.sampling_type['n_iterations']
        else:
            return None

    @property
    def parameter_names(self):
        if isinstance(self.configuration,PyposmatConfigurationFile):
            return [p for p in self.configuration.sampling_distribution]
        else:
            return None

    @property
    def qoi_names(self):
        if isinstance(self.configuration,PyposmatConfigurationFile):
            return [k for k in self.configuration.qois]
        else:
            return None

    @property
    def error_names(self):
        if isinstance(self.configuration,PyposmatConfigurationFile):
            return ['{}.err'.format(k) for k in self.qoi_names]
        else:
            return None

    @property
    def free_parameter_names(self):
        if isinstance(self.configuration,PyposmatConfigurationFile):
            return self.configuration.free_parameter_names

    def initialize_datafile_in(self,filename):
        """ read the datafile

        Args:
            filename(str): path of the filename to read

        """

        assert isinstance(filename,str)

        self.datafile_in_fn = filename

        self.datafile_in = PyposmatDataFile()
        self.datafile_in.read(filename=filename)

    def initialize_datafile_out(self,filename):
        """ configure the datafile for output

        Args:
            filename(str): path of the output file to be written to

        """

        assert isinstance(filename,str)

        self.datafile_out_fn = filename

        self.datafile_out = PyposmatDataFile(filename=filename)

    def initialize_mpi_information(self,mpi_rank=None,mpi_info=None):
        """ initialize mpi information

        Args:
            mpi_rank(int): the rank of the process running this program
            mpi_info(int): the size of the process running this program

        """

        assert isinstance(mpi_rank,int) or mpi_rank is None
        assert isinstance(mpi_size,int) or mpi_size is None

        if mpi_rank is None and mpi_size is None:
            self.mpi_rank = 0
            self.mpi_size = 1 
        else:
            self.mpi_rank = mpi_rank
            self.mpi_size = mpi_size

    def initialize_logfile(self,logfile=None,log_to_stdout=True):
        """ initialize the logger

        Configurtion of the log object has different behavior based upon the type passed
        into the argument o_log.  If o_log is PyposmatLogFile, that object will be accessed
        by reference.  A string is assumed to be a filename location.  By default the
        argument for o_log is None, which means logging will go to standard out by means of 
        the print() function.

        Args:
            logfile (str,PyposmatLogFile,None): default: None
        """

        assert isinstance(log_to_stdout,bool)
        self.log_to_stdout = log_to_stdout

        # if logfile is an instance PyposmatLogFile, the we set the logfile
        # is an attribute by reference
        if type(logfile) is PyposmatLogFile:
            self.logfile = logfile

        # if logfile is a string, then we configure a new instance of the 
        # logfile.
        elif type(logfile) is str:
            self.logfile_fn = logfile
            self.logfile = PyposmatLogFile(filename=self.logfile_fn)

        # the default of setting the logfile
        elif logfile is None:
            self.logfile = None
        else:
            m = "log object must be str, PyposmatLogFile, or None"
            raise TypeError(m)

    def log(self,message):
        """ log a message

        Args:
            message(str): the message to be logged

        """

        if type(message) is str:
            m = message
        elif type(message) is list:
            m = "\n".join(message)
        else:
            raise TypeError()
        
        if isinstance(self.logfile,PyposmatLogFile):
            self.logfile.write(m)

        if self.log_to_stdout:
            print(m)

