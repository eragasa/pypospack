import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from sklearn import manifold
from pypospack.pyposmat.data import datafile


def normalize(df):
    obj_scaler = StandardScaler()
    df = obj_scaler.fit_transform(df)
    return df


def pca_transform(df):
    obj_pca = PCA()
    df = obj_pca.fit_transform(df)
    return df


def tsne_transform(df):
    obj_tsne = manifold.TSNE(init='pca', random_state=None, perplexity=50, n_iter=5000, n_components=3)
    df = obj_tsne.fit_transform(df)
    return df


def cluster(df):
    obj_dbscan = DBSCAN(eps=0.75, min_samples=10)
    labels = obj_dbscan.fit(df).labels_
    return labels


def plot_2d(df, title):
    _clusterids = set(df['cluster_id'])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    print(title)
    for clusterid in _clusterids:
        x = df.loc[df['cluster_id'] == clusterid]['col_0']
        y = df.loc[df['cluster_id'] == clusterid]['col_1']
        print(x.shape, y.shape)
        ax.scatter(x, y, s=1)
    plt.title(s=title)
    plt.show()


def plot_3d(df, title):
    _clusterids = set(df['cluster_id'])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    print(title)
    for clusterid in _clusterids:
        x = df.loc[df['cluster_id'] == clusterid]['col_0']
        y = df.loc[df['cluster_id'] == clusterid]['col_1']
        z = df.loc[df['cluster_id'] == clusterid]['col_2']
        print(x.shape, y.shape, z.shape)
        ax.scatter(x, y, z, s=1)
    plt.title(s=title)
    plt.show()


if __name__ == "__main__":
    datafile = datafile.PyposmatDataFile()
    # datafile.read(r'C:\Users\Seaton\repos\pypospack\visualization_apps\MgO_pareto\data\culled_009.out')
    datafile.read(r'C:\Users\Seaton\repos\pypospack\tests\data_test\Ni__eam__born_exp_fs_00\data__Ni__eam__born_exp_fs_03\pyposmat.kde.10.out')
    # datafile.read(r'C:\Users\Seaton\repos\pypospack\dev\Si_iterative_sampling\data\pyposmat.kde.10.out')
    names = datafile.parameter_names
    data = datafile.df[names]
    nrows, ncols = data.shape
    generic_names = ['col_{}'.format(str(i)) for i in range(ncols)]

    norm_data = normalize(data)
    norm_clusters = cluster(norm_data)
    norm_df = pd.DataFrame(data=norm_data, columns=generic_names)
    norm_df['cluster_id'] = norm_clusters
    plot_3d(norm_df, 'Normalized Only')

    pca_data = pca_transform(norm_data)
    pca_clusters = cluster(pca_data)
    pca_df = pd.DataFrame(data=pca_data, columns=generic_names)
    pca_df['cluster_id'] = pca_clusters
    plot_3d(pca_df, 'Normalized PCA')

    tsne_data = tsne_transform(norm_data)
    tsne_clusters = cluster(tsne_data)
    tsne_df = pd.DataFrame(data=tsne_data, columns=generic_names[:3])
    tsne_df['cluster_id'] = tsne_clusters
    plot_3d(tsne_df, 'Normalized tSNE')
