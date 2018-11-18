from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile


class BasePipeSegment(object):
    """
    Base object for Pyposmat data pipeline objects to inherit from
    """

    def __init__(self, o_logger=None):
        self.o_logger = o_logger  # logging file object
        self.configuration_fn = None
        self.configuration = None
        self.data_fn = None
        self.data = None
        self.df = None

        self.parameter_names = None
        self.error_names = None
        self.qoi_names = None
        self.normalized_names = None  # TODO: break this into param, err, qoi
        self.pca_names = None
        self.manifold_names = None

    def read_configuration(self, filename):
        """
        reads in the pyposmat configuration file
        - set self.configuration
        :param filename: (string)
        """
        self.configuration_fn = filename
        self.configuration = PyposmatConfigurationFile()
        self.configuration.read(filename)

    def read_data(self, filename):
        """
        reads in the pyposmat data filename
        :param filename: (string)
        """
        self.data_fn = filename
        self.data = PyposmatDataFile()
        self.data.read(filename)
        self.df = self.data.df

    def select_data(self, types=['param']):
        """
        selects subset of dataframe
        - mutates self.df
        :param types: list[string]
        :return: pandas.DataFrame
        """
        _names = []
        if 'param' in types:
            _names += self.configuration.parameter_names
        if 'qoi' in types:
            _names += self.configuration.qoi_names
        if 'err' in types:
            _names += self.configuration.error_names
        self.df = self.data.df[_names]
        return self.df

    def log(self, msg):
        if self.o_logger is None:
            print(msg)
        else:
            self.o_logger.write(msg)


class PyposmatPipeline(object):

    def __init__(self):
        pass
