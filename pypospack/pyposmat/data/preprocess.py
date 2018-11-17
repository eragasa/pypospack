from collections import OrderedDict
import pandas as pd
from sklearn import preprocessing
from pypospack.pyposmat.data.pipeline import BasePipeSegment


class PyposmatPreprocessor(BasePipeSegment):

    def __init__(self, o_logger=None):
        super().__init__(o_logger)

    def normalize_data(self, norm_by=None, d=None):
        """
        normalizes the dataframe
        :param d: nested dict of custom args
                  - default to StandardScaler
        :param norm_by: list[string]
                        - df keys to subselect and normalize
                        - normalize entire df if None
        """
        # process arg: norm_by
        if norm_by is None:
            _df = self.df
        else:
            _df = self.df[norm_by]

        # process arg: d
        _kwargs = OrderedDict()
        if d is None:
            self._normalize_standard_scaler(df=_df)  # default normalizer is standard scaler
        else:
            # build kwargs
            for k, v in d['normalizer']['args'].items():
                _kwargs[k] = v
            for k, v in _kwargs.items():
                d['normalizer']['args'][k] = v

            # check normalizer type
            if d['normalizer']['type'] == 'standard_scaler':
                self._normalize_standard_scaler(df=_df, kwargs=_kwargs)  # pass custom args to standard scaler
            else:
                raise NotImplementedError()

    def _normalize_standard_scaler(self, df, kwargs=None):
        """
        mutates self.df by concatenating normalized columns
        """
        # process arg: kwargs
        if kwargs is None:
            normalizer = preprocessing.StandardScaler()
        else:
            normalizer = preprocessing.StandardScaler(**kwargs)

        arr = normalizer.fit_transform(df)
        norm_cols = ["n{}".format(col) for col in list(df)]  # add 'n' to indicate normalization has occurred
        norm_df = pd.DataFrame(data=arr, columns=norm_cols)
        self.df = pd.concat([self.df, norm_df])
