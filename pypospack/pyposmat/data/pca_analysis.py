from collections import OrderedDict
import copy
import pandas as pd

from sklearn.decomposition import PCA

from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile


class PyposmatPcaAnalysis(object):
    def __init__(self, o_logger=None):
        if o_logger is not None:
            self.log = o_logger
        
        self.configuration_fn = None
        self.data_fn = None

        self.data = None
        self.df = None
        self.configuration = None

        self.parameter_names = None
        self.qoi_names = None
        self.error_names = None

    def init_from_ordered_dict(d, o_logger=None):
        """
        constructor from a OrderedDict

        Arguments:
        ==========
        d(collections.OrderedDict):

        """

        o = PyposmatPcaAnalysis(o_logger=o_logger)

        o.configuration_fn = d['configuration_fn']
        try:
            o.read_configuration(filename=o.configuration_fn)
        except FileNotFoundError as e:
            print("ConfigurationError in PyposmatClusterAnalysis")
            for k,v in d.items():
                print("    {}={}".format(k,v))
            print("o.configuration_fn={}".format(
                o.configuration_fn)
            )
            raise e

        o.data_fn = d['data_fn']
        o.read_data(filename=o.data_fn)

        if 'include_parameters' in d:
            o.include_parameters = d['include_parameters']
        else:
            d['include_parameters'] = o.include_parameters
        if 'include_qois' in d:
            o.include_qois = d['include_qois']
        else:
            d['include_qois'] = o.include_qois
        if 'include_errors' in d:
            o.include_errors = d['include_errors']
        else:
            d['include_errors'] = o.include_errors

        if 'normalizer' in d:
            o.normalize_data(**d['normalizer'])

        if 'normalizer_type' in d:
            o.normalizer_type = d['normalizer_type']
            o.normalize_data(o.normalizer_type)
        # quick patch to fix the duplicate axis issue
        if 'cluster_id' in list(o.data.df):
            o.data.df = o.data.df.drop(['cluster_id'], axis=1)
        return o

    # Unresolved 'o' what is o supposed to be?
    # should 'd' be returned?
    def to_dict(self):
        d = OrderedDict()
        d['configuration_fn'] = o.configuration_fn
        d['data_fn'] = o.data_fn
        d['rescaler_type'] = o.rescaler_type
        d['analysis_on']['include_parameters'] = o.include_parameters
        d['analysis_on']['include_errors'] = o.include_errors

    def to_ordered_dict(self):
        d = OrderedDict()
        d['pyposmat_config_fn'] = self.configuration_fn
        d['pyposmat_data_fn'] =self.data_fn

    def read_configuration(self,filename):
        """
        read in pyposmat configuration file

        Argum ents:
        ==========
        filename(str):
        """
        self.configuration_fn = filename
        self.configuration = PyposmatConfigurationFile()
        self.configuration.read(filename=filename)

    def read_data(self,filename):
        """
        read in pyposmat data filename

        Arguments:
        ==========
        filename(str):
        """

        self.data_fn = filename
        self.data = PyposmatDataFile()
        self.data.read(filename)

        self.parameter_names = self.data.parameter_names
        self.qoi_names = self.data.qoi_names
        self.error_names = self.data.error_names
        self.df = self.data.df

    def transform_pca(self, df=None, d=None):
        # process arg: df
        _df = None
        if df is None and self.df is None:
            raise ValueError("no data to transform")
        elif df is None:
            _df = copy.deepcopy(self.df)
        else:
            _df = copy.deepcopy(df)

        # process arg: d
        if d is None:
            return self._transform_pca(_df)  # use default pca with default settings
        else:
            kwargs = OrderedDict()
            for k, v in d['pca']['args'].items():
                kwargs[k] = v
            for k, v in kwargs.items():
                d['pca']['args'][k] = v

            if d['pca']['type'] == 'pca':
                return self._transform_pca(df=_df, kwargs=kwargs)
            elif d['pca']['type'] == 'cca':
                return self._transform_cca(df=_df, kwargs=kwargs)
            elif d['pca']['type'] == 'kernel_pca':
                return self._transform_kernel_pca(df=_df, kwargs=kwargs)

    def _transform_pca(self, df, kwargs=None):
        if kwargs is None:
            o_pca = PCA()
        else:
            o_pca = PCA(**kwargs)

        o_pca.fit(df)
        pca_arr = o_pca.transform(df)
        nrows, ncols = pca_arr.shape
        pca_cols = ["pca_{}".format(i) for i in range(ncols)]
        pca_df = pd.DataFrame(data=pca_arr, columns=pca_cols)
        return pca_df

    def _transform_cca(self, df, kwargs=None):
        raise NotImplementedError()

    def _transform_kernel_pca(self, df, kwargs=None):
        raise NotImplementedError()
