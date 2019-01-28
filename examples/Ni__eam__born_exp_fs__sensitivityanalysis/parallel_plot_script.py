from pypospack.pyposmat.data import PyposmatDataFile, PyposmatConfigurationFile
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN, KMeans
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import pandas as pd
from pandas.tools.plotting import parallel_coordinates

import copy



if __name__ == "__main__":
    configuration_path = "/home/seaton/python-repos/pypospack/examples/Ni__eam__born_exp_fs__sensitivityanalysis/data__from_pareto_optimization/pyposmat.config.in"
    data_path = "/home/seaton/python-repos/pypospack/examples/Ni__eam__born_exp_fs__sensitivityanalysis/data__from_pareto_optimization/pyposmat.kde.5.out"
    
    o_config = PyposmatConfigurationFile()
    o_config.read(configuration_path)

    o_data = PyposmatDataFile()
    o_data.read(data_path)

    o_normalizer = StandardScaler()
    normal_cols = o_normalizer.fit_transform(o_data.parameter_df)
    normal_df = pd.DataFrame(data=normal_cols, columns=o_data.parameter_names)

    o_tsne = TSNE()
    tsne_cols = o_tsne.fit_transform(normal_df)

    # o_cluster = DBSCAN()
    # cluster_ids = o_cluster.fit_predict(tsne_cols)

    o_cluster = KMeans(n_clusters=5)
    cluster_ids = o_cluster.fit_predict(tsne_cols)

    clustered_df = copy.deepcopy(o_data.qoi_df)
    clustered_df['cluster_id'] = cluster_ids

    centroid_df = clustered_df.groupby(['cluster_id']).median()
    centroid_df['cluster_id'] = range(5)

    plt.xticks(rotation=80)
    ax = parallel_coordinates(centroid_df, 'cluster_id', colormap=plt.get_cmap("Set1"))
    plt.title("QOIs Clustered in Parameter Space by tSNE+KMeans")
    plt.tight_layout()
    plt.savefig("parallel_plot_qoi_by_cluster.png")
    plt.show()
