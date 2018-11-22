import pandas as pd
from sklearn import manifold
from collections import OrderedDict
from pypospack.pyposmat.data.pipeline import BasePipeSegment


class PyposmatManifoldAnalysis(BasePipeSegment):

    def __init__(self, o_logger=None):
        super().__init__(o_logger)

    def transform_tsne(self, cols=None, clusters=None, kwargs=None):
        pass

    def calculate_manifold(self, man_by=None, d=None):
        """
        normalizes the dataframe
        :param d: nested dict of custom args
                  - default to TSNE
        :param man_by: list[string]
                        - df keys to learn manifold over
                        - learn entire df if None
        """
        # process arg: man_by
        if man_by is None:
            _df = self.df
        else:
            _df = self.df[man_by]

        # process arg: d
        _kwargs = OrderedDict()
        if d is None:
            self._calculate_manifold_tsne(df=_df)  # default to tsne
        else:
            # build kwargs
            for k, v in d['manifold']['args'].items():
                _kwargs[k] = v
            for k, v in _kwargs.items():
                d['manifold']['args'][k] = v

            # check manifold type
            if d['manifold']['type'] == 'tsne':
                self._calculate_manifold_tsne(df=_df, kwargs=_kwargs)
            else:
                raise NotImplementedError()

    def _calculate_manifold_tsne(self, df, kwargs=None):
        """
        mutates self.df by concatenating tsne columns
        """
        # process arg: kwargs
        if kwargs is None:
            o_tsne = manifold.TSNE()  # default args
        else:
            o_tsne = manifold.TSNE(**kwargs)

        arr = o_tsne.fit_transform(df)
        nrows, ncols = arr.shape
        tsne_cols = ["tsne_{}".format(i) for i in range(ncols)]
        self.manifold_names = tsne_cols
        tsne_df = pd.DataFrame(data=arr, columns=tsne_cols)
        self.df = pd.concat([self.df, tsne_df])
