import pandas as pd
from sklearn import manifold
from collections import OrderedDict
from pypospack.pyposmat.data.pipeline import BasePipeSegment


class PyposmatManifoldAnalysis(BasePipeSegment):

    def __init__(self):
        super().__init__()

    def transform_tsne(self, abs_cols=None, cols=None, clusters=None, kwargs=None):
        # process arg: cols, clusters
        df = self.select_data(cols=cols, clusters=clusters)
        # process arg: abs_cols
        if abs_cols is not None:
            df = df[abs_cols]
        # process arg: kwargs
        kwargs = self.process_kwargs('tsne', kwargs)
        o_tsne = manifold.TSNE(**kwargs)

        arr = o_tsne.fit_transform(df)
        nrows, ncols = arr.shape
        tsne_cols = ['tsne_{}'.format(i) for i in range(ncols)]
        self.manifold_names = tsne_cols
        tsne_df = pd.DataFrame(data=arr, columns=tsne_cols)
        self.df = pd.concat([self.df, tsne_df], axis=1)

    def transform_mds(self, abs_cols=None, cols=None, clusters=None, kwargs=None):
        # process arg: cols, clusters
        df = self.select_data(cols=cols, clusters=clusters)
        # process arg: abs_cols
        if abs_cols is not None:
            df = df[abs_cols]
        # process arg: kwargs
        kwargs = self.process_kwargs('mds', kwargs)
        o_mds = manifold.MDS(**kwargs)

        arr = o_mds.fit_transform(df)
        nrows, ncols = arr.shape
        mds_cols = ['mds_{}'.format(i) for i in range(ncols)]
        self.manifold_names = mds_cols
        mds_df = pd.DataFrame(data=arr, columns=mds_cols)
        self.df = pd.concat([self.df, mds_df], axis=1)

    def transform_spectral(self, abs_cols=None, cols=None, clusters=None, kwargs=None):
        # process arg: cols, clusters
        df = self.select_data(cols=cols, clusters=clusters)
        # process arg: abs_cols
        if abs_cols is not None:
            df = df[abs_cols]
        # process arg: kwargs
        kwargs = self.process_kwargs('spectral', kwargs)
        o_spectral = manifold.SpectralEmbedding(**kwargs)

        arr = o_spectral.fit_transform(df)
        nrows, ncols = arr.shape
        spectral_cols = ['spectral_{}'.format(i) for i in range(ncols)]
        self.manifold_names = spectral_cols
        spectral_df = pd.DataFrame(data=arr, columns=spectral_cols)
        self.df = pd.concat([self.df, spectral_df], axis=1)

    def transform_lle(self, abs_cols=None, cols=None, clusters=None, kwargs=None):
        # process arg: cols, clusters
        df = self.select_data(cols=cols, clusters=clusters)
        # process arg: abs_cols
        if abs_cols is not None:
            df = df[abs_cols]
        # process arg: kwargs
        kwargs = self.process_kwargs('lle', kwargs)
        o_lle = manifold.LocallyLinearEmbedding(**kwargs)

        arr = o_lle.fit_transform(df)
        nrows, ncols = arr.shape
        lle_cols = ['lle_{}'.format(i) for i in range(ncols)]
        self.manifold_names = lle_cols
        lle_df = pd.DataFrame(data=arr, columns=lle_cols)
        self.df = pd.concat([self.df, lle_df], axis=1)

    def transform_isomap(self, abs_cols=None, cols=None, clusters=None, kwargs=None):
        # process arg: cols, clusters
        df = self.select_data(cols=cols, clusters=clusters)
        # process arg: abs_cols
        if abs_cols is not None:
            df = df[abs_cols]
        # process arg: kwargs
        kwargs = self.process_kwargs('isomap', kwargs)
        o_isomap = manifold.Isomap(**kwargs)

        arr = o_isomap.fit_transform(df)
        nrows, ncols = arr.shape
        isomap_cols = ['isomap_{}'.format(i) for i in range(ncols)]
        self.manifold_names = isomap_cols
        isomap_df = pd.DataFrame(data=arr, columns=isomap_cols)
        self.df = pd.concat([self.df, isomap_df], axis=1)
