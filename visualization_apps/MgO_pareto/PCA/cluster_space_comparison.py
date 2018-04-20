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
    obj_tsne = manifold.TSNE(random_state=None, perplexity=50, n_iter=5000, n_components=3)
    df = obj_tsne.fit_transform(df)
    return df


def cluster(df):
    obj_dbscan = DBSCAN(eps=0.75, min_samples=15)
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


def phillpot_plot(df1, df2):
    #df1 is tsne
    #df2 is pca
    _clusterids1 = set(df1['cluster_id'])
    _clusterids2 = set(df2['cluster_id'])
    f = plt.figure('t-SNE Space Visualization of Clusters', figsize=(8,8))
    ax_ul = plt.subplot(221, aspect='equal', adjustable='box-forced')
    ax_ur = plt.subplot(222, aspect='equal', adjustable='box-forced')
    ax_ll = plt.subplot(223, aspect='equal', adjustable='box-forced')
    ax_lr = plt.subplot(224, aspect='equal', adjustable='box-forced')
    '''
    x_ul_min, x12 = ax1.get_xlim()
    y11, y12 = ax1.get_ylim()
    ax1.set_aspect(abs(x12 - x11) / abs(y12 - y11))

    x21, x22 = ax2.get_xlim()
    y21, y22 = ax2.get_ylim()
    ax2.set_aspect(abs(x22 - x11) / abs(y22 - y21))
    '''
    print(_clusterids1)
    print(_clusterids2)
    # tsne with noise
    for clusterid1 in _clusterids1:
        x1 = df1.loc[df1['cluster_id'] == clusterid1]['col_0']
        y1 = df1.loc[df1['cluster_id'] == clusterid1]['col_1']
        ax_ul.scatter(x1, y1, s=1)

    # tsne without noise
    _clusterids1.remove(-1)
    for clusterid1 in _clusterids1:
        x1 = df1.loc[df1['cluster_id'] == clusterid1]['col_0']
        y1 = df1.loc[df1['cluster_id'] == clusterid1]['col_1']
        ax_ll.scatter(x1, y1, s=1)

    # pca with noise
    for clusterid2 in _clusterids2:
        x2 = df2.loc[df2['cluster_id'] == clusterid2]['col_0']
        y2 = df2.loc[df2['cluster_id'] == clusterid2]['col_1']
        ax_ur.scatter(x2, y2, s=1)

    # pca without noise
    _clusterids2.remove(-1)
    for clusterid2 in _clusterids2:
        x2 = df2.loc[df2['cluster_id'] == clusterid2]['col_0']
        y2 = df2.loc[df2['cluster_id'] == clusterid2]['col_1']
        ax_lr.scatter(x2, y2, s=1)

    ax_ul.set_title('t-SNE Clustered')
    ax_ur.set_title('PCA Clustered')
    plt.show()


if __name__ == "__main__":
    datafile = datafile.PyposmatDataFile()
    # datafile.read(r'C:\Users\Seaton\repos\pypospack\visualization_apps\MgO_pareto\data\culled_009.out')
    datafile.read(r'/home/seaton/repos/pypospack/tests/data_test/Ni__eam__born_exp_fs_00/data__Ni__eam__born_exp_fs_03/pyposmat.kde.10.out')
    # datafile.read(r'C:\Users\Seaton\repos\pypospack\dev\Si_iterative_sampling\data\pyposmat.kde.10.out')
    names = datafile.parameter_names
    data = datafile.df[names]
    nrows, ncols = data.shape
    generic_names = ['col_{}'.format(str(i)) for i in range(ncols)]

    norm_data = normalize(data)
    norm_clusters = cluster(norm_data)
    # transform to tsne space
    norm_tsne_arr = tsne_transform(norm_data)
    norm_tsne_df = pd.DataFrame(data=norm_tsne_arr, columns=generic_names[:3])
    norm_tsne_df['cluster_id'] = norm_clusters
    #plot_3d(norm_tsne_df, 'Normalized Only in t-SNE Space')

    pca_data = pca_transform(norm_data)
    pca_clusters = cluster(pca_data)
    # transform to tsne space
    pca_tsne_arr = tsne_transform(pca_data)
    pca_tsne_df = pd.DataFrame(data=pca_tsne_arr, columns=generic_names[:3])
    pca_tsne_df['cluster_id'] = pca_clusters
    #plot_3d(pca_tsne_df, 'Normalized PCA in t-SNE Space')

    tsne_data = tsne_transform(norm_data)
    tsne_clusters = cluster(tsne_data)
    tsne_df = pd.DataFrame(data=tsne_data, columns=generic_names[:3])
    tsne_df['cluster_id'] = tsne_clusters
    #plot_3d(tsne_df, 'Normalized tSNE')

    phillpot_plot(tsne_df, pca_tsne_df)
