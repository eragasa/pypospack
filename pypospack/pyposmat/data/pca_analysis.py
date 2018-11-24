from collections import OrderedDict
import pandas as pd
from sklearn.decomposition import PCA, FastICA, CCA, KernelPCA
from pypospack.pyposmat.data.pipeline import BasePipeSegment


class PyposmatPcaAnalysis(BasePipeSegment):
    def __init__(self):
        super().__init__()

    def transform_pca(self, abs_cols=None, cols=None, clusters=None, kwargs=None):
        # process arg: cols, clusters
        df = self.select_data(cols=cols, clusters=clusters)
        # process arg: abs_cols
        if abs_cols is not None:
            df = df[abs_cols]
        # process arg: kwargs
        kwargs = self.process_kwargs('pca', kwargs)
        o_pca = PCA(**kwargs)

        arr = o_pca.fit_transform(df)
        nrows, ncols = arr.shape
        pca_cols = ["pca_{}".format(i) for i in range(ncols)]
        self.pca_names = pca_cols
        pca_df = pd.DataFrame(data=arr, columns=pca_cols)
        self.df = pd.concat([self.df, pca_df], axis=1)

    def transform_ica(self, abs_cols=None, cols=None, clusters=None, kwargs=None):
        # process arg: cols, clusters
        df = self.select_data(cols=cols, clusters=clusters)
        # process arg: abs_cols
        if abs_cols is not None:
            df = df[abs_cols]
        # process arg: kwargs
        kwargs = self.process_kwargs('ica', kwargs)
        o_ica = FastICA(**kwargs)

        arr = o_ica.fit_transform(df)
        nrows, ncols = arr.shape
        ica_cols = ["ica_{}".format(i) for i in range(ncols)]
        self.pca_names = pca_cols
        ica_df = pd.DataFrame(data=arr, columns=ica_cols)
        self.df = pd.concat([self.df, ica_df], axis=1)

    def transform_cca(self, abs_cols=None, cols=None, clusters=None, kwargs=None):
        # process arg: cols, clusters
        df = self.select_data(cols=cols, clusters=clusters)
        # process arg: abs_cols
        if abs_cols is not None:
            df = df[abs_cols]
        # process arg: kwargs
        kwargs = self.process_kwargs('cca', kwargs)
        o_cca = CCA(**kwargs)

        arr = o_cca.fit_transform(df)
        nrows, ncols = arr.shape
        cca_cols = ["cca_{}".format(i) for i in range(ncols)]
        self.pca_names = cca_cols
        cca_df = pd.DataFrame(data=arr, columns=cca_cols)
        self.df = pd.concat([self.df, cca_df], axis=1)

    def transform_kernel_pca(self, abs_cols=None, cols=None, clusters=None, kwargs=None):
        # process arg: cols, clusters
        df = self.select_data(cols=cols, clusters=clusters)
        # process arg: abs_cols
        if abs_cols is not None:
            df = df[abs_cols]
        # process arg: kwargs
        kwargs = self.process_kwargs('kernel_pca', kwargs)
        o_kernel_pca = KernelPCA(**kwargs)

        arr = o_kernel_pca.fit_transform(df)
        nrows, ncols = arr.shape
        kernel_pca_cols = ["kernel_pca_{}".format(i) for i in range(ncols)]
        self.pca_names = kernel_pca_cols
        kernel_pca_df = pd.DataFrame(data=arr, columns=kernel_pca_cols)
        self.df = pd.concat([self.df, kernel_pca_df], axis=1)
