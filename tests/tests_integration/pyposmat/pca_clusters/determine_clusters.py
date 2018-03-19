import numpy as np
import pandas as pd

from sklearn import preprocessing
from sklearn import metrics
from sklearn import cluster
from sklearn import decomposition
from sklearn.datasets.samples_generator import make_blobs
from pypospack.pyposmat.data import PyposmatDataFile

class PcaClusterer(object):

    def __init__(self):
        self.configuration_fn = None
        self.configuration = None

        self.datafile_fn = None
        self.datafile = None

        self.preprocessor = None
        self.decomposer = None
        self.clusterer = None

    def read_configuration_file(self,filename):
        self.configuration = PyposmatConfigurationFile()
        self.configuration.read(filename=filename)
        self.configuration_fn = filename
    
    def read_datafile(self,filename):
        self.datafile = PyposmatDatafile()
        self.datafile.read(filename=filename)
        self.datafile_fn = filename

    def run(self):
        pass
class ParetoSurface(object):

    def __init__(self):
        self.config = None
        self.data = None:

if __name__ == "__main__":



    # #############################################################################
    # Generate sample data
    centers = [[1, 1], [-1, -1], [1, -1]]
    X, labels_true = make_blobs(n_samples=750, centers=centers, cluster_std=0.4,
                                random_state=0)

    X = preprocessing.StandardScaler().fit_transform(X)
    print(labels_true)
    print(X.shape)
    # #############################################################################
    # Compute DBSCAN
    db = cluster.DBSCAN(eps=0.3, min_samples=10).fit(X)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

    print('Estimated number of clusters: %d' % n_clusters_)
    print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
    print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
    print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))
    print("Adjusted Rand Index: %0.3f"
          % metrics.adjusted_rand_score(labels_true, labels))
    print("Adjusted Mutual Information: %0.3f"
          % metrics.adjusted_mutual_info_score(labels_true, labels))
    print("Silhouette Coefficient: %0.3f"
          % metrics.silhouette_score(X, labels))

    fn = "../../../data_test/Ni__eam__born_exp_bjs_00/data__Ni__eam__born_exp_bjs_01/pyposmat.results.0.mod.out"
    print("loading {}...".format(fn))
    data = PyposmatDataFile()
    data.read(filename=fn)
    nr,nc = data.df.shape
    print("there are {} rows".format(nr))
    
    names = data.parameter_names + data.error_names
    X = data.df[names]
    nr,nc = X.shape
    preprocessor=preprocessing.StandardScaler()
    X_n = preprocessor.fit_transform(X)
    
    pca = decomposition.PCA(n_components=3)
    pca.fit(X_n)
    X_pca = pca.transform(X_n)
    print(pca.explained_variance_ratio_)
    #print(pca.singular_values_)

    print(names[0])
    print(names.remove(names[0]))
    dbscan_eps=0.3
    min_samples=10
    clusterer=cluster.DBSCAN(eps=dbscan_eps,min_samples=10).fit(X_pca)
    core_samples_mask = np.zeros_like(clusterer.labels_, dtype=bool)
    core_samples_mask[clusterer.core_sample_indices_] = True
    labels = clusterer.labels_
    
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    print('Estimated number of clusters: %d' % n_clusters_)
    #print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
    #print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
    #print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))
    #print("Adjusted Rand Index: %0.3f"
    #  % metrics.adjusted_rand_score(labels_true, labels))
    #print("Adjusted Mutual Information: %0.3f"
    #  % metrics.adjusted_mutual_info_score(labels_true, labels))
    print("Silhouette Coefficient: %0.3f" % metrics.silhouette_score(X_pca, labels))



    print(preprocessor.mean(axis=0))
    print(preprocessor.std(axis=0))

    

    exit()
    # #############################################################################
    # Plot result
    import matplotlib.pyplot as plt

    # Black removed and is used for noise instead.
    unique_labels = set(labels)
    colors = [plt.cm.Spectral(each)
              for each in np.linspace(0, 1, len(unique_labels))]
    for k, col in zip(unique_labels, colors):
        if k == -1:
            # Black used for noise.
            col = [0, 0, 0, 1]

        class_member_mask = (labels == k)

        xy = X[class_member_mask & core_samples_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
                 markeredgecolor='k', markersize=14)

        xy = X[class_member_mask & ~core_samples_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
                 markeredgecolor='k', markersize=6)

    plt.title('Estimated number of clusters: %d' % n_clusters_)
    #plt.show()
