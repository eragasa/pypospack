import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering, AffinityPropagation
from sklearn.preprocessing import StandardScaler, normalize
import copy

from pypospack.pyposmat.data import datafile


class TransformPCA(object):
    '''
    0) use StandardScaler to normalize the data by removing mean and scaling to unit variance
    1) transform the df to PCA space
    2) group across the entire df with the {DBSCAN} method
    3) plot the clusters
    '''

    def __init__(self, df, col_names, normalizer, cluster_by, fit_by):
        self.df = df
        self.col_names = col_names
        self.normalizer = normalizer
        self.cluster_by = cluster_by
        self.fit_by = fit_by

        self.norm_pca_df = self.normalize_and_transform()
        self.norm_pca_df['cluster_id'] = self.make_clusters()

    def normalize_and_transform(self):
        if self.fit_by == 'all':
            df = copy.deepcopy(self.df)
        elif self.fit_by == 'param':
            df = self.df[self.col_names['param']]
        elif self.fit_by == 'err':
            df = self.df[self.col_names['err']]
        elif self.fit_by == 'qoi':
            df = self.df[self.col_names['qoi']]
        else:
            df = None
            print('unsupported fit by method')
            exit()

        print(df)

        if self.normalizer == 'sphere':
            df = normalize(df)
        elif self.normalizer == 'standard':
            obj_normalize = StandardScaler()
            df = obj_normalize.fit_transform(df)
        else:
            df = None
            print("unsupported normalize method")
            exit()
        nrows, ncols = df.shape
        # this naming convention is wrong and i feel bad about doing it
        names = ['pca_{}'.format(i) for i in range(ncols)]
        transformed_df = pd.DataFrame(data=df, columns=names)
        return transformed_df

    def make_clusters(self):
        if self.cluster_by == 'kmeans':
            obj_kmeans = KMeans(n_clusters=3)
            obj_kmeans = obj_kmeans.fit(self.norm_pca_df)
            labels = obj_kmeans.labels_
        elif self.cluster_by == 'dbscan':
            obj_dbscan = DBSCAN()
            obj_dbscan = obj_dbscan.fit(self.norm_pca_df)
            labels = obj_dbscan.labels_
        else:
            labels = None
            print("unsupported clustering method")
            exit()
        return labels

    def plot_3d(self):
        _clusterids = set(self.norm_pca_df['cluster_id'])
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for clusterid in _clusterids:
            x = self.norm_pca_df.loc[self.norm_pca_df['cluster_id'] == clusterid]['pca_0']
            y = self.norm_pca_df.loc[self.norm_pca_df['cluster_id'] == clusterid]['pca_1']
            z = self.norm_pca_df.loc[self.norm_pca_df['cluster_id'] == clusterid]['pca_2']
            print(x.shape, y.shape, z.shape)
            ax.scatter(x, y, z, s=1)
        plt.title(s='PCA Transformation')
        plt.show()


if __name__ == "__main__":
    datafile = datafile.PyposmatDataFile()
    datafile.read(r'/home/seaton/repos/pypospack/visualization_apps/MgO_pareto/data/culled_009.out')
    data = datafile.df
    names = {'err': datafile.error_names, 'param': datafile.parameter_names, 'qoi': datafile.qoi_names}
    t_pca = TransformPCA(df=data, col_names=names, normalizer='sphere', cluster_by='kmeans', fit_by='all')
    t_pca.plot_3d()
