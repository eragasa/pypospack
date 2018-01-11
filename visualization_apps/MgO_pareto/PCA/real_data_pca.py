import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import copy

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
    error_df = df[error_names]
    qoi_df = df[qoi_names]

    return {'total_df': df, 'param_df': parameter_df, 'err_df': error_df, 'qoi_df': qoi_df}, \
           {'param_names': parameter_names, 'err_names': error_names, 'qoi_names': qoi_names}

def cluster(results):
    # results is tuple of dicts from read_data_file
    data_dict = results[0]
    total_df = data_dict['total_df']
    total_df = total_df.drop('sim_id', axis=1)
    # using arbitrary 5 clusters for now
    kmeans = KMeans(n_clusters=5)
    # determine the clusters of each space and adds cluster index for later grouping
    total_df['cluster_index'] = kmeans.fit_predict(total_df)

    return total_df


def plot_pca(results):
    total_df = cluster(results)
    # first group the total by cluster_index
    groupings = []
    # 5 because 5 clusters
    for i in range(5):
        group = total_df.loc[total_df['cluster_index'] == i]
        groupings.append(group)
    group_indices = {str(i): groups.index.tolist() for i, groups in enumerate(groupings)}
    # group_indices structured {'0': [index1, index2...], '1': [index1...], ...}

    # do pca on each spatial grouping (param, err, qoi)
    # plot results
    colors = ['red', 'orange', 'yellow', 'green', 'blue']
    pca = PCA()
    pca_info_dict = {}
    data_dict = results[0]
    param_pca_array = pca.fit_transform(data_dict['param_df'])
    pca_info_dict['param'] = {'explained_variance': pca.explained_variance_ratio_,
                              'vectors': param_pca_array}
    err_pca_array = pca.fit_transform(data_dict['err_df'])
    pca_info_dict['err'] = {'explained_variance': pca.explained_variance_ratio_,
                            'vectors': err_pca_array}
    qoi_pca_array = pca.fit_transform(data_dict['qoi_df'])
    pca_info_dict['err'] = {'explained_variance': pca.explained_variance_ratio_,
                            'vectors': qoi_pca_array}
    print(pca_info_dict['param']['explained_variance'])

    for array, name in zip([param_pca_array, err_pca_array, qoi_pca_array],\
                        ['Parameter Space', 'Error Space', 'QOI Space']):
        plt.figure()
        for cluster_id, c in zip(group_indices, colors):
            relevant_indices = group_indices[cluster_id]
            plt.scatter(array[relevant_indices, 0],
                        array[relevant_indices, 1],
                        color=c, label=cluster_id)
        plt.legend(loc='best', shadow=False, scatterpoints=1)
        plt.title('PCA of '+name)
    plt.show()


if __name__ == "__main__":
    results = read_data_file(r'C:\Users\Seaton\repos\pypospack\visualization_apps\MgO_pareto\data\culled_009.out')
    plot_pca(results)


