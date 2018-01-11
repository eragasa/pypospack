import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

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


if __name__ == "__main__":
    results = read_data_file(r'C:\Users\Seaton\repos\pypospack\visualization_apps\MgO_pareto\data\culled_009.out')

    # analyze paramaters
    param_df = results[0]['param_df']
    from sklearn.preprocessing import normalize
    param_df = normalize(param_df)
    obj_pca = PCA()
    pca_param_np = obj_pca.fit_transform(param_df)
    nrows, ncols = pca_param_np.shape
    pca_param_names = ['param_pca_{}'.format(i) for i in range(ncols)]
    pca_param_df = pd.DataFrame(data=pca_param_np, columns=pca_param_names)

    # create df of vectors and values to write excel file easily
    eigen_values = obj_pca.explained_variance_
    eigen_vectors = obj_pca.components_
    eigen_df = pd.DataFrame(columns=pca_param_names)
    eigen_vectors = eigen_vectors.tolist()
    for i, names in enumerate(pca_param_names):
        eigen_df[names] = eigen_vectors[i]
        if i == len(eigen_vectors) - 1:
            eigen_df.loc[i + 1] = eigen_values.tolist()
    writer = pd.ExcelWriter('eigen_table.xlsx')
    EIGEN_TABLE = eigen_df.to_excel(writer, 'Sheet 1')
    writer.save()

    # cluster the PCA transformed data
    obj_kmeans = KMeans(n_clusters=5)
    obj_kmeans = obj_kmeans.fit(pca_param_df)
    labels = obj_kmeans.labels_
    pca_param_df['cluster_id'] = labels

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    _clusterids = set(pca_param_df['cluster_id'])
    # fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    for clusterid in _clusterids:
        x = pca_param_df.loc[pca_param_df['cluster_id'] == clusterid]['param_pca_1']
        y = pca_param_df.loc[pca_param_df['cluster_id'] == clusterid]['param_pca_2']
        z = pca_param_df.loc[pca_param_df['cluster_id'] == clusterid]['param_pca_3']
        plt.scatter(x,y,s=1)
    plt.show()
    