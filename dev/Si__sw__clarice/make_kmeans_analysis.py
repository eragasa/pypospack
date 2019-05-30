"""
example on how to GMM

adapted from
Jake VanderPlass, Python Data Science Handbook
https://jakevdp.github.io/PythonDataScienceHandbook/05.12-gaussian-mixtures.html
"""
import os
import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import numpy as np
import pandas as pd

# generate data
from sklearn.datasets.samples_generator import make_blobs
X, y_true = make_blobs(n_samples=400,centers=4,cluster_std=0.60,random_state=0)
X = X[:,::-1]
print('type(X):{}'.format(type(X)))
# generate clusters
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
n_kmeans_clusters = 4
kmeans = KMeans(n_kmeans_clusters, random_state=0)
print(type(kmeans))
cluster_id = kmeans.fit(X).predict(X)
cluster_ids= set(cluster_id)
#print(cluster_id)
plt.scatter(X[:, 0], X[:, 1], c=cluster_id, s=40, cmap='viridis')
plt.show()

def plot_kmeans(kmeans,X,n_clusters=4,r_seed=0,ax=None):
    assert isinstance(kmeans,KMeans)
    assert isinstance(X,np.ndarray)
    assert isinstance(n_clusters,int)
    assert isinstance(r_seed,int)
    assert ax is None \
            or isinstance(ax,plt.Axes)

    cluster_id = kmeans.fit_predict(X)

    if ax is None:
        fig, ax = plt.subplots(1,1)
    
    ax.scatter(X[:,0],X[:,1],c=cluster_id, s=40, cmap='viridis')

    centers = kmeans.cluster_centers_
    radii = [cdist(X[cluster_id == i], [center]).max()
                 for i, center in enumerate(centers)]

    for c, r in zip(centers, radii):
        ax.add_patch(plt.Circle(c, r, fc='#CCCCCC', lw=3, alpha=0.5, zorder=1))

    plt.show()

kmeans = KMeans(n_clusters=4, random_state=0)
plot_kmeans(kmeans, X)
