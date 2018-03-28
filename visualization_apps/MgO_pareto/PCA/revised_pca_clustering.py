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

    def __init__(self, df, normalizer, cluster_by):
        self.df = df
        self.normalizer = normalizer
        self.cluster_by = cluster_by
        self.pca_df = self.normalize_and_transform()

    def normalize_and_transform(self):
        df = copy.deepcopy(self.df)
        if self.normalizer == 'sphere':
            df = normalize(df)
        elif self.normalizer == 'standard':
            obj_normalize = StandardScaler()
            df = obj_normalize.fit_transform(df)
        else:
            df = None
            print("unsupported normalize method")
            exit()
        obj_pca = PCA()
        pca_df = obj_pca.fit_transform(df)
        self.write_pca_data(obj_pca)
        nrows, ncols = pca_df.shape
        # this naming convention is wrong and i feel bad about doing it
        names = ['pca_{}'.format(i) for i in range(ncols)]
        pca_df = pd.DataFrame(data=pca_df, columns=names)
        labels = self.make_clusters(pca_df)
        pca_df['cluster_id'] = labels
        return pca_df

    def make_clusters(self, pca_df):
        if self.cluster_by == 'kmeans':
            obj_kmeans = KMeans(n_clusters=3)
            obj_kmeans = obj_kmeans.fit(pca_df)
            labels = obj_kmeans.labels_
        elif self.cluster_by == 'dbscan':
            obj_dbscan = DBSCAN()
            obj_dbscan = obj_dbscan.fit(pca_df)
            labels = obj_dbscan.labels_
        else:
            labels = None
            print("unsupported clustering method")
            exit()
        return labels

    def write_pca_data(self, obj_pca):
        eigen_vectors = obj_pca.components_
        eigen_values = obj_pca.explained_variance_
        vector_df = pd.DataFrame(data=eigen_vectors, columns=['pca_{}'.format(str(i)) for i in range(eigen_vectors.shape[1])])
        value_df = pd.DataFrame(data=eigen_values, columns=['eigen_values'])
        result = pd.concat([value_df, vector_df], axis=1)
        result.to_excel('param_err_sphere_dbscan_pca.xlsx')

    def plot_3d(self):
        _clusterids = set(self.pca_df['cluster_id'])
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for clusterid in _clusterids:
            x = self.pca_df.loc[self.pca_df['cluster_id'] == clusterid]['pca_0']
            y = self.pca_df.loc[self.pca_df['cluster_id'] == clusterid]['pca_1']
            z = self.pca_df.loc[self.pca_df['cluster_id'] == clusterid]['pca_2']
            print(x.shape, y.shape, z.shape)
            ax.scatter(x, y, z, s=1)
        plt.title(s='Normalizer: {n} Cluster: {c} '.format(n=self.normalizer, c=self.cluster_by))
        plt.show()


if __name__ == "__main__":
    datafile = datafile.PyposmatDataFile()
    datafile.read(r'/home/seaton/repos/pypospack/visualization_apps/MgO_pareto/data/culled_009.out')
    param_data = datafile.parameter_df
    err_data = datafile.error_df
    data = pd.concat([param_data, err_data], axis=1)
    t_pca = TransformPCA(df=data, normalizer='sphere', cluster_by='dbscan')
    t_pca.plot_3d()
