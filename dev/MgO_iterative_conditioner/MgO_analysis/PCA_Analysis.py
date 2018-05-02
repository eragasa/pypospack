import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pypospack.pyposmat.data import PyposmatDataFile
from sklearn import datasets
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

#load data
data_fn = "pyposmat.results.0.out"
data = PyposmatDataFile()
data.read(filename=data_fn)


# print(data.df['MgO_NaCl.p11'])
data.df = pd.concat(
    [
        data.df,
        pd.DataFrame(
            data.df['MgO_NaCl.p11'].abs(),
            columns = ['MgO_NaCl.p11.abs']
        )
    ],
    axis = 0
)
# data.df = data.df.dropna(axis=0)
print(
    data.df
)


# print(data.df.nsmallest(30, 'MgO_NaCl.p11.abs'))

pca = PCA(n_components=2)
pca.fit(data.df[data.parameter_names])

X_50000 = pca.transform(
   data.df.nsmallest(50000, 'MgO_NaCl.p11')
)
X_25000 = pca.transform(
   data.df.nsmallest(25000, 'MgO_NaCl.p11')
)
X_10000 = pca.transform(
   data.df.nsmallest(10000, 'MgO_NaCl.p11')
)
X_10000 = pca.transform(
   data.df.nsmallest(10000, 'MgO_NaCl.p11')
)



plt.figure()
# colors = ['navy', 'turquoise', 'darkorange']
lw = 2
plt.scatter(X_50000[:,0], X_50000[:,1],label='50000')
plt.scatter(X_25000[:,0], X_25000[:,1],label='25000')
plt.scatter(X_10000[:,0], X_10000[:,1],label='10000')
plt.title('PCA of dataset')


plt.show()
