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
