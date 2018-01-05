import matplotlib.pyplot as plt
from sklearn import datasets
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

iris = datasets.load_iris()
colors = ['navy', 'turquoise', 'darkorange']
lw = 2

X = iris.data
y = iris.target
target_names = iris.target_names

pca = PCA(n_components=2)
X_r = pca.fit(X).transform(X)

plt.figure()
# add kmeans to pre-pca plot
kmeans = KMeans(n_clusters=3)
kmeans = kmeans.fit(X)
centroids = kmeans.cluster_centers_
for color, i, target_name in zip(colors, [0, 1, 2], target_names):
    plt.scatter(X[y == i, 0], X[y == i, 1], color=color, alpha=.8, lw=lw,
                label=target_name)
    plt.scatter(centroids[i][0], centroids[i][1], color='black')
plt.legend(loc='best', shadow=False, scatterpoints=1)
plt.title('Pre-PCA of IRIS dataset')


# Percentage of variance explained for each components
print('explained variance ratio (first two components): %s'
      % str(pca.explained_variance_ratio_))

plt.figure()
# add kmeans to pca plot
kmeans = KMeans(n_clusters=3)
kmeans = kmeans.fit(X_r)
centroids = kmeans.cluster_centers_
for color, i, target_name in zip(colors, [0, 1, 2], target_names):
    plt.scatter(X_r[y == i, 0], X_r[y == i, 1], color=color, alpha=.8, lw=lw,
                label=target_name)
    plt.scatter(centroids[i][0], centroids[i][1], color='black')
plt.legend(loc='best', shadow=False, scatterpoints=1)
plt.title('PCA of IRIS dataset')

plt.show()
