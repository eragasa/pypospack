import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, DBSCAN
from sklearn.manifold import TSNE
from pypospack.pyposmat.data.datafile import PyposmatDataFile
from pypospack.kde import chiu1999_h, silverman1986_h


def angle_between(a1, a2):
    a1_u = a1/np.linalg.norm(a1)
    a2_u = a2/np.linalg.norm(a2)
    return np.arccos(np.clip(np.dot(a1_u, a2_u), -1.0, 1.0))


# load the kde file
filename = "/home/seaton/repos/pypospack/dev/cluster_sampling/dev__cluster_sampling/data/pyposmat.kde.0.out"
data = PyposmatDataFile(filename=filename)
data.read(filename)

# normalize the data
scaler = StandardScaler()
df = data.parameter_df
df = scaler.fit_transform(df)
df = pd.DataFrame(data=df, columns=data.parameter_names)

# learn the tSNE manifold of normalized parameters
tsne = TSNE()
tsne_dims = tsne.fit_transform(df)

# cluster the tsne space
# dbscan = DBSCAN()
# labels = dbscan.fit_predict(tsne_dims)
kmeans = KMeans(n_clusters=10)
kmeans.fit(tsne_dims)
labels = kmeans.labels_

# concat the labels with the df
df['cluster_id'] = labels

# split the df by cluster
cluster_ids = set(labels)
subselections = []
for c in cluster_ids:
    sub = df.loc[df['cluster_id'] == c]
    # append only the parameters
    sub = sub[data.parameter_names]
    subselections.append(sub)

# check bandwidth of each cluster
for s in subselections:
    print("Chiu Bandwidth: {c}\nSilverman Bandwidth: {s}\n".format(c=chiu1999_h(s.values.T),
                                                                   s=silverman1986_h(s.values.T)))

print("\n")

# calculate the angle difference between the primary principle component for each cluster compared to the first
standard_vector = None
for s in subselections:
    pca = PCA()
    pca.fit(s)
    first_component = pca.components_[0]
    if standard_vector is None:
        standard_vector = first_component
        print("Angle Difference: 0")
    else:
        angle = angle_between(standard_vector, first_component)
        print("Angle Difference: {}".format(angle))
