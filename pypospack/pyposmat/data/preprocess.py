from collections import OrderedDict
import pandas as pd
from sklearn import preprocessing
from pypospack.pyposmat.data import BasePipeSegment


class PyposmatPreprocessor(BasePipeSegment):

    def __init__(self):
        super().__init__()

    def normalize_standard_scaler(self, abs_cols=None, cols=None, clusters=None, kwargs=None):
        # process arg: cols, clusters
        df = self.select_data(cols=cols, clusters=clusters)
        # process arg: abs_cols
        if abs_cols is not None:
            df = df[abs_cols]  
        # process arg: kwargs
        kwargs = self.process_kwargs('standard_scaler', kwargs)
        o_normalizer = preprocessing.StandardScaler(**kwargs)

        arr = o_normalizer.fit_transform(df)
        norm_cols = ["n{}".format(col) for col in list(df)]
        norm_df = pd.DataFrame(data=arr, columns=norm_cols)
        self.df = pd.concat([self.df, norm_df], axis=1)  # mutate self.df
        self._set_names(cols=cols, clusters=clusters)  # mutate normalized names

    def _set_names(self, cols, clusters):
        if 'param' in cols:
            self.n_parameter_names = ["n{}".format(name) for name in self.parameter_names]
        if 'err' in cols:
            self.n_error_names = ["n{}".format(name) for name in self.error_names]
        if 'qoi' in cols:
            self.n_qoi_names = ["n{}".format(name) for name in self.qoi_names]
