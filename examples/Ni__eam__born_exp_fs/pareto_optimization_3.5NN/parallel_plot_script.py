from pypospack.pyposmat.data import PyposmatDataFile, PyposmatConfigurationFile
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN, KMeans
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import pandas as pd
from pandas.tools.plotting import parallel_coordinates

import copy



if __name__ == "__main__":
    # define paths to configuration and data files
    configuration_path = "/home/seaton/python-repos/pypospack/examples/Ni__eam__born_exp_fs/pareto_optimization_3.5NN/data/new_data/pyposmat.config.in"
    data_path = "/home/seaton/python-repos/pypospack/examples/Ni__eam__born_exp_fs/pareto_optimization_3.5NN/data/new_data/pyposmat.kde.4.out"

    # init the configuration object
    o_config = PyposmatConfigurationFile()
    o_config.read(configuration_path)

    # init the data file object
    o_data = PyposmatDataFile()
    o_data.read(data_path)

    # normalize the QOIs to prepare for tSNE
    o_normalizer = StandardScaler()
    normal_qoi_arr = o_normalizer.fit_transform(o_data.qoi_df)
    normal_qoi_df = pd.DataFrame(data=normal_qoi_arr, columns=o_data.qoi_names)

    # learn and apply tSNE manifold to the normal QOIs
    o_tsne = TSNE()
    tsne_arr = o_tsne.fit_transform(normal_qoi_df)

    # find KMeans clusters in the tSNE space
    o_cluster = KMeans(n_clusters=3)
    cluster_ids = o_cluster.fit_predict(tsne_arr)

    # plot the tSNE space by cluster
    normal_qoi_df['tsne_0'] = tsne_arr[:,0]
    normal_qoi_df['tsne_1'] = tsne_arr[:,1]
    normal_qoi_df['cluster_id'] = cluster_ids
    for cid in set(cluster_ids):
        subset = normal_qoi_df.loc[normal_qoi_df['cluster_id'] == cid]
        _x = subset['tsne_0']
        _y = subset['tsne_1']
        plt.scatter(x=_x, y=_y, s=1, cmap=plt.get_cmap("Set1"))
    plt.savefig("qoi_clusters_in_tsne_space.png")
    plt.show()

    # generate a DataFrame containing percent errors
    percent_error_df = copy.deepcopy(o_data.qoi_df)
    for q_name in o_data.qoi_names:
        percent_error_df[q_name] = (o_data.qoi_df[q_name]-o_config.qoi_targets[q_name])/o_config.qoi_targets[q_name]
    percent_error_df['cluster_id'] = cluster_ids

    print("min: ", o_data.qoi_df["Ni_fcc.isf"].min())
    print("max: ", o_data.qoi_df["Ni_fcc.isf"].max())
    print("mean: ", o_data.qoi_df["Ni_fcc.isf"].mean())
    print("target: ", o_config.qoi_targets["Ni_fcc.isf"])
    # centroid_df = percent_error_df.groupby(['cluster_id']).median()
    # centroid_df['cluster_id'] = range(5)

    plt.xticks(rotation=80)
    ax = parallel_coordinates(percent_error_df, 'cluster_id', colormap=plt.get_cmap("Set1"))
    plt.title("Percent Error Across QOI Clusters (tSNE+KMeans)")
    plt.tight_layout()
    plt.savefig("parallel_plot_pct_err_by_qoi_cluster.png")
    plt.show()
