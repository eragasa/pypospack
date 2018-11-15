from collections import OrderedDict

import pandas as pd

from sklearn import preprocessing

from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile


class PyposmatPreprocessor(object):

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

    def read_configuration(self, filename):
        """
        read in pyposmat configuration file

        Argum ents:
        ==========
        filename(str):
        """
        self.configuration_fn = filename
        self.configuration = PyposmatConfigurationFile()
        self.configuration.read(filename=filename)

    def read_data(self, filename):
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

    def select_data(self, types=['param']):
        _names = []
        if 'param' in types:
            _names += self.configuration.parameter_names
        if 'qoi' in types:
            _names += self.configuration.qoi_names
        if 'err' in types:
            _names += self.configuration.error_names

        self.df = self.data.df[_names]
        return self.df

    def normalize_data(self, df=None, d=None):
        # process arg: df
        _df = None
        if df is None and self.df is None:
            raise ValueError("no data to normalize")
        elif df is None:
            _df = self.df
        else:
            _df = df

        # process arg: d
        if d is None:
            return self._normalize_standard_scaler(df=_df)
        elif d['normalizer']['type'] == 'standard_scaler':
            return self._normalize_standard_scaler(df=_df, d=d)
        else:
            raise NotImplementedError()

    def _normalize_standard_scaler(self, df, d=None):
        # process arg: d
        if d is None:
            normalizer = preprocessing.StandardScaler()
        else:
            kwargs = OrderedDict()
            for k, v in d['normalizer']['args'].items():
                kwargs[k] = v

            for k, v in kwargs.items():
                d['normalizer']['args'][k] = v

            normalizer = preprocessing.StandardScaler(**kwargs)

        arr = normalizer.fit_transform(df)
        norm_cols = ["n{}".format(col) for col in list(df)]  # add 'n' to indicate normalization
        norm_df = pd.DataFrame(data=arr, columns=norm_cols)
        return norm_df
