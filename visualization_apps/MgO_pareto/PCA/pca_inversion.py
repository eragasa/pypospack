import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize
from sklearn.cluster import KMeans, DBSCAN

from pypospack.pyposmat.data import datafile


def make_clusters(pca_df, cluster_by):
    if cluster_by == 'kmeans':
        obj_kmeans = KMeans(n_clusters=10)
        obj_kmeans = obj_kmeans.fit(pca_df)
        labels = obj_kmeans.labels_
    elif cluster_by == 'dbscan':
        obj_dbscan = DBSCAN(eps=0.75, min_samples=10).fit(pca_df)
        labels = obj_dbscan.labels_
    else:
        labels = None
        print("unsupported clustering method")
        exit()
    return labels


def plot_3d(df):
    _clusterids = set(df['cluster_id'])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for clusterid in _clusterids:
        x = df.loc[df['cluster_id'] == clusterid]['inv_0']
        y = df.loc[df['cluster_id'] == clusterid]['inv_1']
        z = df.loc[df['cluster_id'] == clusterid]['inv_2']
        print(x.shape, y.shape, z.shape)
        ax.scatter(x, y, z, s=1)
    plt.title(s='Inverted Space Cluster Projection')
    plt.show()


def compare_inversion(df):
    # start with raw data df and do transform and inversion in here
    # write raw df to csv
    df.to_csv(path_or_buf='raw_data.csv', sep=',')
    orig_np, orig_norms = normalize(df, return_norm=True)
    nrows, ncols = orig_np.shape
    orig_names = ['orig_{}'.format(i) for i in range(ncols)]
    orig_df = pd.DataFrame(data=orig_np, columns=orig_names)

    obj_pca = PCA()
    pca_np = obj_pca.fit_transform(orig_df)
    pca_names = ['pca_{}'.format(i) for i in range(ncols)]
    pca_df = pd.DataFrame(data=pca_np, columns=pca_names)

    labels = make_clusters(pca_df, 'kmeans')

    # undo pca to get to normalized data
    inv_np = obj_pca.inverse_transform(pca_np)
    inv_names = ['inv_{}'.format(i) for i in range(ncols)]
    inv_df = pd.DataFrame(data=inv_np, columns=inv_names)

    assert orig_df.shape == pca_df.shape == inv_df.shape

    # undo normalization to get to raw data
    for i, n in enumerate(orig_norms.tolist()):
        inv_df.iloc[i] = inv_df.iloc[i] * n

    # compare error between raw input and inverted output
    difference = np.array(inv_df) - np.array(df)

    norm = stats.norm(scale=difference.std(), loc=difference.mean())
    x = np.linspace(difference.min(), difference.max(), 100)
    plt.plot(x, norm.pdf(x), 'r-', lw=5, alpha=0.6, label='pdf of errors')
    plt.hist(difference)
    plt.legend(loc='best')
    plt.title(s='Error Analysis of PCA Inversion')
    plt.xlabel(s='error')
    plt.ylabel(s='frequency')
    plt.show()

    # include the cluster labels
    inv_df['cluster_id'] = labels
    # write denormalized df to csv
    inv_df.to_csv(path_or_buf='fully_inverted_data.csv', sep=',')

    plot_3d(inv_df)




if __name__ == "__main__":
    datafile = datafile.PyposmatDataFile()
    datafile.read(r'/home/seaton/repos/pypospack/visualization_apps/MgO_pareto/data/culled_009.out')
    #datafile.read('/home/seaton/repos/pypospack/tests/data_test/Ni__eam__born_exp_fs_00/data__Ni__eam__born_exp_fs_03/pyposmat.kde.10.out')
    #datafile.read('/home/seaton/repos/pypospack/dev/Si_iterative_sampling/data/pyposmat.kde.10.out')
    names = datafile.parameter_names
    data = datafile.df[names]
    compare_inversion(data)
