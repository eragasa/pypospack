import copy
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data.pipeline_configuration import PyposmatPipelineConfigurationFile


class BasePipeSegment(object):
    """
    Base object for Pyposmat data pipeline objects to inherit from
    """

    def __init__(self):
        self.o_logger = None
        self.configuration = None
        self.df = None

        self.parameter_names = None
        self.error_names = None
        self.qoi_names = None
        self.n_parameter_names = None  # normalized
        self.n_error_names = None  # normalized
        self.n_qoi_names = None  # normalized
        self.pca_names = None
        self.manifold_names = None

    def process_kwargs(key, d):
        if d is None:
            return {}  # use default args
        try:
            kwargs = d[key]  # use configred args
        except KeyError:
            kwargs = {}  # use default args
        return kwargs

    def select_data(self, cols=None, clusters=None):
        # no subselection if both are none
        if cols is None and clusters is None:
            return self.df

        _df = copy.deepcopy(self.df)
        if clusters is not None:  # subselect by cluster ids
            df_list = []
            for cid in clusters:
                _df = self.df.loc[self.df['cluster_id'] == cid]
                df_list.append(_df)
            _df = pd.concat(df_list)
        if cols is not None:  # subselect by column types
            _names = []
            if 'param' in cols:
                _names += self.parameter_names
            if 'qoi' in cols:
                _names += self.qoi_names
            if 'err' in cols:
                _names += self.error_names
            if 'n_param' in cols:
                _names += self.n_parameter_names
            if 'n_qoi' in cols:
                _names += self.n_qoi_names
            if 'n_err' in cols:
                _names += self.n_error_names
            if 'pca' in cols:
                _names += self.pca_names
            if 'manifold' in cols:
                _names += self.manifold_names
            _df = _df[_names]
        return _df

    def log(self, msg):
        if self.o_logger is None:
            print(msg)
        else:
            self.o_logger.write(msg)


class PyposmatPipeline(object):

    def __init__(self,
                 o_logger=None,
                 configuration_fn=None,
                 data_fn=None,
                 df=None):
        self.o_logger = o_logger  # logging file object
        self.configuration = self.read_configuration(configuration_fn)
        self.data = None
        if df is not None:
            self.df = df
            self.parameter_names = None
            self.error_names = None
            self.qoi_names = None
        elif data_fn is not None:
            self.read_data(data_fn)
            self.df = data.df
            self.parameter_names = data.parameter_names
            self.error_names = data.error_names
            self.qoi_names = data.qoi_names
        else:
            raise ValueError("no data to work with")
            exit()

        self.n_parameter_names = None  # normalized
        self.n_error_names = None  # normalized
        self.n_qoi_names = None  # normalized
        self.pca_names = None
        self.manifold_names = None

    def read_configuration(self, filename):
        pass

    def read_data(self, filename):
        self.data = PyposmatDataFile()
        self.data.read(filename)

    def reset(self, segment):
        assert type(segment) is BasePipeSegment
        self.df = segment.df
        self.parameter_names = segment.parameter_names
        self.error_names = segment.error_names
        self.qoi_names = segment.qoi_names
        self.n_parameter_names = segment.n_parameter_names
        self.n_error_names = segment.n_error_names
        self.n_qoi_names = segment.n_qoi_names
        self.pca_names = segment.pca_names
        self.manifold_names = segment.manifold_names
