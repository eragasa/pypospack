import copy
import yaml
import time
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.io.filesystem import OrderedDictYAMLLoader


class BasePipeSegment(object):
    """
    Base object for Pyposmat data pipeline objects to inherit from
    """

    def __init__(self):
        self.o_logger = None
        self.df = None

        self.parameter_names = None
        self.error_names = None
        self.qoi_names = None
        self.n_parameter_names = None  # normalized
        self.n_error_names = None  # normalized
        self.n_qoi_names = None  # normalized
        self.pca_names = None
        self.manifold_names = None

    def process_kwargs(self, key, d):
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
        self.configuration_fn = configuration_fn
        self.configuration = None
        self.data_fn = data_fn
        self.data = None
        self.df = df

        self.parameter_names = None
        self.error_names = None
        self.qoi_names = None
        self.n_parameter_names = None  # normalized
        self.n_error_names = None  # normalized
        self.n_qoi_names = None  # normalized
        self.pca_names = None
        self.manifold_names = None

    def read_configuration(self, filename):
        with open(filename, 'r') as f:
            config = yaml.load(f, OrderedDictYAMLLoader)
        self.configuration = config

    def write_configuration(self, filename, d):
        with open(filename, 'w') as f:
            yaml.dump(d, f, default_flow_style=False)

    def read_data(self, filename):
        self.data = PyposmatDataFile()
        self.data.read(filename)
        self.df = self.data.df
        self.parameter_names = self.data.parameter_names
        self.error_names = self.data.error_names
        self.qoi_names = self.data.qoi_names

    def log(self, msg):
        if self.o_logger is None:
            print(msg)
        else:
            self.o_logger.write(msg)

    def reset(self, o_segment):
        assert isinstance(o_segment, BasePipeSegment)
        self.df = o_segment.df
        self.parameter_names = o_segment.parameter_names
        self.error_names = o_segment.error_names
        self.qoi_names = o_segment.qoi_names
        self.n_parameter_names = o_segment.n_parameter_names
        self.n_error_names = o_segment.n_error_names
        self.n_qoi_names = o_segment.n_qoi_names
        self.pca_names = o_segment.pca_names
        self.manifold_names = o_segment.manifold_names

    def spawn_pipeline_segment(self, segment_type):
        if segment_type == 'preprocess':
            from pypospack.pyposmat.data.preprocess import PyposmatPreprocessor
            o_segment = PyposmatPreprocessor()
        elif segment_type == 'pca':
            from pypospack.pyposmat.data.pca_analysis import PyposmatPcaAnalysis
            o_segment = PyposmatPcaAnalysis()
        elif segment_type == 'cluster':
            from pypospack.pyposmat.data.cluster_analysis import SeatonClusterAnalysis
            o_segment = SeatonClusterAnalysis()
        elif segment_type == 'manifold':
            from pypospack.pyposmat.data.manifold_analysis import PyposmatManifoldAnalysis
            o_segment = PyposmatManifoldAnalysis()
        elif segment_type == 'plot':
            from pypospack.pyposmat.data.plotting import PyposmatPlotter
            o_segment = PyposmatPlotter()
        else:
            raise ValueError("unknown segment type")

        o_segment.o_logger = self.o_logger
        o_segment.df = self.df
        o_segment.parameter_names = self.parameter_names
        o_segment.error_names = self.error_names
        o_segment.qoi_names = self.qoi_names
        o_segment.n_parameter_names = self.n_parameter_names
        o_segment.n_error_names = self.n_error_names
        o_segment.n_qoi_names = self.n_qoi_names
        o_segment.pca_names = self.pca_names
        o_segment.manifold_names = self.manifold_names

        return o_segment

    def make_function_calls(self, o_segment, calls):
        for index in calls:
            self.log("calling function {}".format(calls[index]['function']))
            func = getattr(o_segment, calls[index]['function'])
            kwargs = calls[index]['args']
            func(**kwargs)

    def run(self):
        pipeline_start_time = time.time()
        for index in self.configuration:
            self.log("starting step {} of {}".format(index+1, len(self.configuration)))  # +1 to count like a normal person
            step_start_time = time.time()
            o_segment = self.spawn_pipeline_segment(self.configuration[index]['segment_type'])
            self.make_function_calls(o_segment=o_segment,
                                     calls=self.configuration[index]['function_calls'])
            self.reset(o_segment)
            step_end_time = time.time()
            step_delta = step_end_time-step_start_time
            step_delta = round(step_delta, 4)
            self.log("step {} complete in {} seconds".format(index+1, step_delta))
        pipeline_end_time = time.time()
        pipeline_delta = pipeline_end_time-pipeline_start_time
        pipeline_delta = round(pipeline_delta, 4)
        self.log("pipeline complete in {} seconds\n".format(pipeline_delta))
