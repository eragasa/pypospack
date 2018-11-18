from collections import OrderedDict
import pandas as pd
from sklearn.decomposition import PCA
from pypospack.pyposmat.data.pipeline import BasePipeSegment


class PyposmatPcaAnalysis(BasePipeSegment):
    def __init__(self, o_logger=None):
        super().__init__(o_logger)

    def transform_pca(self, pca_by=None, d=None):
        """
        transforms dataframe to pca space
        :param pca_by: list[string]
                       - df keys to subselect and transform
                       - transform entire df is None
        :param d: nested dict of custom args
                  - default to sklearn PCA implementation
        """
        # process arg: pca_by
        if pca_by is None:
            _df = self.df
        else:
            _df = self.df[pca_by]

        # process arg: d
        _kwargs = OrderedDict()
        if d is None:
            self._transform_pca(df=_df)  # default pca implementation
        else:
            # build kwargs
            for k, v in d['pca']['args'].items():
                _kwargs[k] = v
            for k, v in _kwargs.items():
                d['pca']['args'][k] = v

            # check pca type
            if d['pca']['type'] == 'pca':
                self._transform_pca(df=_df, kwargs=_kwargs)  # custom args to PCA
            elif d['pca']['type'] == 'cca':
                self._transform_cca(df=_df, kwargs=_kwargs)  # custom args to CCA
            elif d['pca']['type'] == 'kernel_pca':
                self._transform_kernel_pca(df=_df, kwargs=_kwargs)  # custom args to Kernel PCA
            else:
                raise NotImplementedError()

    def _transform_pca(self, df, kwargs=None):
        """
        mutates self.df by concatentating pca columns
        """
        # process arg: kwargs
        if kwargs is None:
            o_pca = PCA()
        else:
            o_pca = PCA(**kwargs)

        arr = o_pca.fit_transform(df)
        nrows, ncols = arr.shape
        pca_cols = ["pca_{}".format(i) for i in range(ncols)]
        self.pca_names = pca_cols
        pca_df = pd.DataFrame(data=arr, columns=pca_cols)
        self.df = pd.concat([self.df, pca_df])

    def _transform_cca(self, df, kwargs=None):
        raise NotImplementedError()

    def _transform_kernel_pca(self, df, kwargs=None):
        raise NotImplementedError()
