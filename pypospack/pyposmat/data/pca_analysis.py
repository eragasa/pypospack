from collections import OrderedDict
import pandas as pd
from sklearn.decomposition import PCA
from pypospack.pyposmat.data.pipeline import BasePipeSegment


class PyposmatPcaAnalysis(BasePipeSegment):
    def __init__(self):
        super().__init__()

    def transform_pca(self, cols=None, clusters=None, kwargs=None):
        # process arg: cols, clusters
        df = self.select_data(cols=cols, clusters=clusters)
        # process arg: kwargs
        kwargs = self.process_kwargs('pca', kwargs)
        o_pca = PCA(**kwargs)

        arr = o_pca.fit_transform(df)
        nrows, ncols = arr.shape
        pca_cols = ["pca_{}".format(i) for i in range(ncols)]
        self.pca_names = pca_cols
        pca_df = pd.DataFrame(data=arr, columns=pca_cols)
        self.df = pd.concat([self.df, pca_df])

    def transform_cca(self, cols=None, clusters=None, kwargs=None):
        raise NotImplementedError()

    def transform_kernel_pca(self, cols=None, clusters=None, kwargs=None):
        raise NotImplementedError()
