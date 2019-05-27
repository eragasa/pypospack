import numpy as np
import pandas as pd
import scipy.stats
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import normalize
import matplotlib.pyplot as plt


def read_data_file(filename):

    with open(filename, 'r') as f:
        lines = f.readlines()

    names = [s.strip() for s in lines[0].strip().split(',')]
    types = [s.strip() for s in lines[1].strip().split(',')]

    all_values = []
    for i in range(2, len(lines)):
        line = lines[i].strip()
        values = [float(s.strip()) for s in line.split(',')]
        values[0] = int(values[0])
        all_values.append(list(values))
    values = np.array(all_values)

    parameter_names = [
        n for i, n in enumerate(names) \
        if types[i] == 'param']
    qoi_names = [
        n for i, n in enumerate(names) \
        if types[i] == 'qoi']
    error_names = [
        n for i, n in enumerate(names) \
        if types[i] == 'err']
    df = pd.DataFrame(data=values, columns=names, copy=True)
    df.set_index('sim_id')
    parameter_df = df[parameter_names]
    # static params are: chrg_O, MgMg_A, MgMg_rho, MgMg_C, MgO_C
    parameter_df = parameter_df.drop(['chrg_O', 'MgMg_A', 'MgMg_rho', 'MgMg_C', 'MgO_C'], axis=1)
    # make rho exponential to put on equal standing with other factors
    parameter_df['MgO_rho'] = np.exp(parameter_df['MgO_rho'])
    parameter_df['OO_rho'] = np.exp(parameter_df['OO_rho'])
    error_df = df[error_names]
    qoi_df = df[qoi_names]

    return {'total_df': df, 'param_df': parameter_df, 'err_df': error_df, 'qoi_df': qoi_df}, \
           {'param_names': parameter_names, 'err_names': error_names, 'qoi_names': qoi_names}


def sample_by_pca(n_samples, params_in, n_clusters):
    '''cluster data by PCA, split clusters, do kde sampling in each cluster, return samples in real space'''

    # normalize and store norms
    df, norms = normalize(params_in, return_norm=True, axis=0)

    # convert to PCA space
    obj_pca = PCA()
    pca_np = obj_pca.fit_transform(df)
    nrows, ncols = pca_np.shape
    pca_names = ['pca_{}'.format(i) for i in range(ncols)]
    pca_df = pd.DataFrame(data=pca_np, columns=pca_names)

    # do kmeans clustering
    obj_kmeans = KMeans(n_clusters=n_clusters)
    obj_kmeans = obj_kmeans.fit(pca_df)
    labels = obj_kmeans.labels_

    # add cluster labels column
    pca_df['cluster_id'] = labels

    # split groups into individual pandas frames
    cluster_df_list = [pca_df.loc[pca_df['cluster_id'] == cluster] for cluster in set(labels)]

    # do the resampling in normalized pca space
    resampled_clusters = {}
    for cluster_id, clusters in enumerate(cluster_df_list):
        clusters = clusters.drop(['cluster_id'], axis=1)
        resampled_cols = []
        for cols in clusters.columns:
            kde_sampler = scipy.stats.gaussian_kde(dataset=clusters[cols])
            new_samples = kde_sampler.resample(n_samples)
            new_col = pd.Series(data=new_samples[0], name=cols+'_resample')
            resampled_cols.append(new_col)
        resample_pca_df = pd.concat(resampled_cols, axis=1)

        # invert to non pca space
        resample_inv_np = obj_pca.inverse_transform(resample_pca_df)
        resample_inv_names = ['resample_inv_{}'.format(i) for i in range(len(resampled_cols))]
        resample_inv_df = pd.DataFrame(data=resample_inv_np, columns=resample_inv_names)

        # undo original normalization
        for c, n in zip(resample_inv_df.columns, norms):
            resample_inv_df.loc[:, c] *= n
        resampled_clusters['cluster_' + str(cluster_id)] = resample_inv_df

    return resampled_clusters


if __name__ == "__main__":
    data_file = r'/home/seaton/repos/pypospack/visualization_apps/MgO_pareto/data/culled_009.out'
    data, names = read_data_file(filename=data_file)
    # for proper de-normalization len samples needs to be len input
    # possibly change order of operations to solve
    new_samples = sample_by_pca(n_samples=len(data['param_df']), params_in=data['param_df'], n_clusters=3)

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.hist(data['param_df']['chrg_Mg'], bins=50)
    ax1.set_title('original distribution of chrg_Mg')
    ax2 = fig.add_subplot(212)
    ax2.hist(new_samples['cluster_0']['resample_inv_0'], bins=50)
    ax2.set_title('resampled distribution of chrg_Mg')
    fig.tight_layout(h_pad=1)
    plt.show()
