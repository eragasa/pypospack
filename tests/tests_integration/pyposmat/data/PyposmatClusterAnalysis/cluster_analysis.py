from collections import OrderedDict
from sklearn import manifold

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import Imputer


from pypospack.data import PyposmatDataFile
from pypospack.data import PyposmatConfigurationFile

class PyposmatClusterSampler(object):

    def __init__(self):
        pass

class PyposmatPreprocessingPipeline(object):
    def __init__(self):
        pass


class PyposmatClusterAnalysis(object):
    def __init__(self):
        self.configuration_fn = None
        self.data_fn = None

        self._data = None
        self._configuration = None
        self.pre_pipeline = None
        self. = None

    def init_from_json(self,json):
        pass

    def init_from_ordered_dict(d):
        """
        constructor from a OrderedDict

        Arguments:
        ==========
        d(collections.OrderedDict):

        """

        o = PyposmatClusterAnalysis()

        o.configuration_fn = d['configuration_fn']
        o.read_configuration(filename=self.configuration_fn)

        o.data_fn = d['pyposmat_data_fn']
        o.read_datafile(filename=self.data_fn)

        if 'normalizer' in d:
            o.normalize_data(**d['normalizer'])

        if 'normalizer_type' in d:
            o.normalizer_type = d['normalizer_type']
            o.normalize_data(o.normalizer_type)

        return o

    def normalize_data(normalizer_type,**kwargs):d
        if normalizer

    def to_dict(self):
        d = OrderedDict()
        d['configuration_fn'] = o.configuration_fn
        d['data_fn'] = o.data_fn
        d['normalizer_type'] = o.




    def to_json(self):
        pass

    def to_ordered_dict(self):
        d = OrderedDict()
        d['pyposmat_config_fn'] = self.configuration_fn
        d['pyposmat_data_fn'] =self.data_fn


    @property
    def configuration(self):
        return self._configuration

    @configuation.setter
    def configuration(self.config):
        self._configuration = config

    def read_configuration(self,filename):
        """
        read in pyposmat configuration file

        Arguments:
        ==========
        filename(str):
        """
        self.configuration_fn = filename
        self.configuration = PyposmatConfigurationFile()
        self.configuration.read(filename=filename)

    def read_datafile(self,filename):
        """
        read in pyposmat data filename

        Arguments:
        ==========
        filename(str):
        """

        self.datafile_fn = filename
        self.datafile = PyposmatDataFile()
        self.datafile.read(filename)

        self.parameter_names = self.datafile.parameter_names
        self.qoi_names = self.datafile.paramter_names
        self.error_names = self.datafile.error_names
