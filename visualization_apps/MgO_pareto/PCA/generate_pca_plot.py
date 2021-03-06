import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

def fake_data(num_groups, shape):
    assert len(shape) == 2
    # create mean and covariance matrix for each group
    means = [create_mean(shape[1]) for i in range(num_groups)]
    covs = [create_cov(shape[1]) for i in range(num_groups)]
    y = create_target(num_groups, shape[0])

    distributions = []
    for m, c in zip(means, covs):
        dist = np.random.multivariate_normal(m, c, shape[0]//num_groups)
        distributions.append(dist)
    distributions = tuple(distributions)
    final_array = np.concatenate(distributions)

    names = ['group'+str(i) for i in range(num_groups)]

    return {'data': final_array, 'target': y, 'names': names}



def create_target(num_groups, samples):
    y = [i for i in range(num_groups) for r in range(samples//num_groups)]
    y = np.array(y)
    return y

def create_mean(dimension):
    mean = [np.random.rand() for i in range(dimension)]
    return mean

def create_cov(dimension):
    cov = []
    for i in range(dimension):
        sub_cov = [0 for i in range(dimension)]
        sub_cov[i] = np.random.rand()
        cov.append(sub_cov)
    return cov

def raw_to_pca(data_array, target_array, target_names):
    assert len(data_array.shape) == 2
    assert len(target_array) == data_array.shape[0]

    pca = PCA()
    # the error on transform does not seem to be an issue
    pca_array = pca.fit(data_array).transform(data_array)
    print('explained variance ratio: ' + str(pca.explained_variance_ratio_))

    # create a second data array which is just a transformed version of the first
    # use to compare with non transformed separation
    transformed_array = data_array*2
    transformed_pca_array = pca.fit(transformed_array).transform(transformed_array)
    print('explained variance ratio (transformed): ' + str(pca.explained_variance_ratio_))

    colors = ['red', 'blue', 'purple', 'black', 'brown', 'yellow', 'orange', 'green']
    colors = colors[:len(target_names)]
    kmeans = KMeans(n_clusters=len(target_names))

    # create figure of the non transformed array
    plt.figure()
    kmeans = kmeans.fit(pca_array)
    centroids = kmeans.cluster_centers_
    # plotting the first 2 dimensions
    for i, c, n in zip(range(len(target_names)), colors, target_names):
        plt.scatter(pca_array[target_array == i, 0], pca_array[target_array == i, 1], alpha=.8, color=c, label=n)
        plt.scatter(centroids[i][0], centroids[i][1], color='black', label=n)
    plt.legend(loc='best', shadow=False, scatterpoints=1)
    plt.title('PCA Test')

    # create figure of the transformed array
    plt.figure()
    kmeans = kmeans.fit(transformed_pca_array)
    centroids = kmeans.cluster_centers_
    # plotting the first 2 dimensions
    for i, c, n in zip(range(len(target_names)), colors, target_names):
        plt.scatter(transformed_pca_array[target_array == i, 0], transformed_pca_array[target_array == i, 1], alpha=.8, color=c, label=n)
        plt.scatter(centroids[i][0], centroids[i][1], color='black', label=n)
    plt.legend(loc='best', shadow=False, scatterpoints=1)
    plt.title('Transformed PCA Test')

    plt.show()

    # look for distributions to be the same in transformed and non transformed plots

if __name__ == "__main__":
    # later pass in actual data from sims but now test with general data
    test_data = fake_data(2, (100, 3))
    raw_to_pca(test_data['data'], test_data['target'], test_data['names'])