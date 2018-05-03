from collections import OrderedDict
import pandas as pd
import numpy as np

from scipy import stats
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

# imports for graphics
import matplotlib.pyplot as plt

from pypospack.pyposmat.data import PyposmatDataFile

#load data
data_fn = "pyposmat.results.0.out"
data = PyposmatDataFile()
data.read(filename=data_fn)
data.df = pd.concat(
    [
        data.df,
        pd.DataFrame(
            data.df['MgO_NaCl.p11'].abs().as_matrix(),
            columns = ['MgO_NaCl.p11.abs']
        ),
        pd.DataFrame(
            data.df['MgO_NaCl.a0.err'].abs().as_matrix(),
            columns =['MgO_NaCl.a0.nerr']
        )
    ],
    axis = 1
)

pca =PCA(n_components=2)
pca.fit(data.df[data.parameter_names])
_pca = pca.transform(data.df[data.parameter_names])
_pca_1_min = _pca[:,0].min()
_pca_1_max = _pca[:,0].max()
_pca_2_min = _pca[:,1].min()
_pca_2_max = _pca[:,1].max()
sample_sizes = [50000,10000,5000,1000]

pca_p11_data=OrderedDict()
for i in sample_sizes:
    df = data.df.nsmallest(i, 'MgO_NaCl.p11.abs')
    pca_p11_data[i] = pca.transform(df[data.parameter_names])

pca_a0_data=OrderedDict()
for i in sample_sizes:
    df = data.df.nsmallest(i, 'MgO_NaCl.a0.nerr')
    pca_a0_data[i] = pca.transform(df[data.parameter_names])

# colors = ['navy', 'turquoise', 'darkorange']
n_plots = len(sample_sizes)
fig, ax = plt.subplots(n_plots,2)
i = 0
#X,Y = np.mgrid[
#    _pca_1_min:_pca_1_max:100j,
#    _pca_2_min:_pca_2_max:100j
#]
#XY = np.vstack([X.ravel(),Y.ravel()])
for k,v in pca_p11_data.items():
    #kde = stats.gaussian_kde(pca_p11_data[k].T)
    #kde_Z = np.reshape(
    #    kde(XY).T,
    #    X.shape
    #)
    #ax[i][0].imshow(
    #    np.rot90(kde_Z),
    #    extent = [_pca_1_min,_pca_1_max,_pca_2_min,_pca_2_max]
    #)
    ax[i][0].scatter(
        pca_p11_data[k][:,0],
        pca_p11_data[k][:,1],
        s=.1)
    i+=1

i = 0
for k,v in pca_a0_data.items():
    #kde = stats.gaussian_kde(pca_a0_data[k].T)
    #kde_Z = np.reshape(
    #    kde(XY).T,
    #    X.shape
    #)
    #ax[i][1].imshow(
    #    np.rot90(kde_Z),
    #    extent = [_pca_1_min,_pca_1_max,_pca_2_min,_pca_2_max]
    #)

    ax[i][1].scatter(
        pca_a0_data[k][:,0],
        pca_a0_data[k][:,1],
        s=.1)
    i+=1

plt.show()

print(data.parameter_names)
print(pca.components_)
print(pca.explained_variance_)
print(pca.explained_variance_ratio_)
