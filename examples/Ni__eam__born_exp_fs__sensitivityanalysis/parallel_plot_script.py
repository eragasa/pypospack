from pypospack.pyposmat.data import PyposmatDataFile, PyposmatConfigurationFile
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN, KMeans
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import pandas as pd
from pandas.tools.plotting import parallel_coordinates

import copy



if __name__ == "__main__":
    configuration_path = "/Users/seaton/python-repos/pypospack/examples/Ni__eam__born_exp_fs__sensitivityanalysis/data__from_pareto_optimization/pyposmat.config.in"
    data_path = "/Users/seaton/python-repos/pypospack/examples/Ni__eam__born_exp_fs__sensitivityanalysis/data__from_pareto_optimization/pyposmat.kde.5.out"

    o_config = PyposmatConfigurationFile()
    o_config.read(configuration_path)

    o_data = PyposmatDataFile()
    o_data.read(data_path)

    o_normalizer = StandardScaler()
    normal_parameter_arr = o_normalizer.fit_transform(o_data.parameter_df)
    normal_parameter_df = pd.DataFrame(data=normal_parameter_arr, columns=o_data.parameter_names)

    o_tsne = TSNE()
    tsne_arr = o_tsne.fit_transform(normal_parameter_df)

    # o_cluster = DBSCAN()
    # cluster_ids = o_cluster.fit_predict(tsne_cols)

    o_cluster = KMeans(n_clusters=5)
    cluster_ids = o_cluster.fit_predict(tsne_arr)

    normal_parameter_df['tsne_0'] = tsne_arr[:,0]
    normal_parameter_df['tsne_1'] = tsne_arr[:,1]
    normal_parameter_df['cluster_id'] = cluster_ids
    for cid in set(cluster_ids):
        subset = normal_parameter_df.loc[normal_parameter_df['cluster_id'] == cid]
        _x = subset['tsne_0']
        _y = subset['tsne_1']
        plt.scatter(x=_x, y=_y, s=1, cmap=plt.get_cmap("Set1"))
    plt.savefig("parameter_clusters_in_tsne_space.png")
    plt.show()


    percent_error_df = copy.deepcopy(o_data.error_df)
    for e_name in o_data.error_names:
        percent_error_df[e_name] = o_data.error_df[e_name].div(o_data.error_df[e_name].mean())
    percent_error_df['cluster_id'] = cluster_ids


    centroid_df = percent_error_df.groupby(['cluster_id']).median()
    centroid_df['cluster_id'] = range(5)

    plt.xticks(rotation=80)
    ax = parallel_coordinates(centroid_df, 'cluster_id', colormap=plt.get_cmap("Set1"))
    plt.title("Normalized Errors Clustered in Parameter Space by tSNE+KMeans")
    plt.tight_layout()
    plt.savefig("parallel_plot_norm_err_by_cluster.png")
    plt.show()
