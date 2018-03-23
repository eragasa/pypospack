import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering, AffinityPropagation
from sklearn.preprocessing import StandardScaler

from pypospack.pyposmat.data import datafile


class TransformPCA(object):
    '''
    0) use StandardScaler to normalize the data by removing mean and scaling to unit variance
    1) transform the df to PCA space
    2) group across the entire df with the {DBSCAN} method
    3) plot the clusters
    '''

    def __init__(self, df):
        self.df = df
        self.norm_pca_df = self.normalize_and_transform()
        self.norm_pca_df['cluster_id'] = self.cluster_by_kmeans()

    def normalize_and_transform(self):
        obj_normalize = StandardScaler()
        df = obj_normalize.fit_transform(self.df)
        nrows, ncols = df.shape
        # this naming convention is wrong and i feel bad about doing it
        names = ['pca_{}'.format(i) for i in range(ncols)]
        transformed_df = pd.DataFrame(data=df, columns=names)
        return transformed_df

    def cluster_by_dbscan(self):
        obj_dbscan = DBSCAN()
        obj_dbscan = obj_dbscan.fit(self.norm_pca_df)
        labels = obj_dbscan.labels_
        return labels

    def cluster_by_kmeans(self):
        obj_kmeans = KMeans(n_clusters=3)
        obj_kmeans = obj_kmeans.fit(self.norm_pca_df)
        labels = obj_kmeans.labels_
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
    datafile.read(r'C:\Users\Seaton\repos\pypospack\visualization_apps\MgO_pareto\data\culled_009.out')
    df = datafile.df
    t_pca = TransformPCA(df)
    t_pca.plot_3d()

