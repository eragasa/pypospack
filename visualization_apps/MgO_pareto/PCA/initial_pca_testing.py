import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA


y = [i for i in [0, 1] for r in range(500)]
y = np.array(y)
names = ['d1', 'd2']

mean_1 = (1, 1, 1)
cov_1 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
x_1 = np.random.multivariate_normal(mean_1, cov_1, 500)

mean_2 = (3, 3, 3)
cov_2 = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
x_2 = np.random.multivariate_normal(mean_2, cov_2, 500)

# build the samples independently to control variance and combine here
x = np.concatenate((x_1, x_2))

plt.figure()
colors = ['red', 'blue']
for i, c, n in zip([0, 1], colors, names):
    plt.scatter(x[y == i, 0], x[y == i, 1], alpha=.8, color=c, label=n)

plt.legend(loc='best', shadow=False, scatterpoints=1)
plt.title('Pre-PCA Test')


pca = PCA()
x_r = pca.fit(x).transform(x)

# Percentage of variance explained for each components
print('explained variance ratio: '+str(pca.explained_variance_ratio_))

plt.figure()
colors = ['red', 'blue']
for i, c, n in zip([0, 1], colors, names):
    plt.scatter(x_r[y == i, 0], x_r[y == i, 1], alpha=.8, color=c, label=n)

plt.legend(loc='best', shadow=False, scatterpoints=1)
plt.title('PCA Test')

plt.show()
