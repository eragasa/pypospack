import warnings

import os
import numpy as np
import numpy.linalg as linalg
import pandas as pd

from scipy import stats

from sklearn import preprocessing
from sklearn.exceptions import DataConversionWarning
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=DataConversionWarning)
import pypospack.utils

from pypospack.pyposmat.data import PyposmatDataFile,PyposmatConfigurationFile

class GaussianKde(gaussian_kde):

    def __init__(self, X, bw_method=None, weights=None):
        self.scaler = preprocessing.StandardScaler().fit(X.T)
        stats.gaussian_kde.__init__(self,self.scaler.transform(X.T).T)

    def evaluate(self, X):
        X_ = self.scaler.transform(np.array(X).T)
        return stats.gaussian_kde.evaluate(self,X_)
     
    def resample(self, size=None):
        X = stats.gaussian_kde.resample(self, size)
        X = self.scaler.inverse_transform(X.T)
        return X.T

if __name__ == "__main__":
    pypospack_path = pypospack.utils.get_pypospack_root_directory()
    config_fn = os.path.join(
            pypospack_path,
            'data','Si__sw__data','pareto_optimization_unconstrained',
            'pyposmat.config.in')
    data_fn = os.path.join(
            pypospack_path,
            'data','Si__sw__data','pareto_optimization_unconstrained',
            'pyposmat.results.1.out')

    configuration = PyposmatConfigurationFile()
    configuration.read(filename=config_fn)

    data = PyposmatDataFile()
    data.read(filename=data_fn)

    qoi_data_df = data.df[configuration.qoi_names]
    n_samples, n_variables = qoi_data_df.shape

    print('n_samples:{}'.format(n_samples))
    print('n_vars:{}'.format(n_variables))
    cov = np.cov(qoi_data_df.T)
    assert cov.shape == (n_variables,n_variables)
    
    eig_values, eig_vectors = linalg.eig(cov)
    eig = zip(eig_values, eig_vectors)

    for k in eig:
        print(k)

    # rescale data
    X = qoi_data_df
    scaler = preprocessing.StandardScaler().fit(X)
    X_scaled = scaler.transform(X)

    cov = np.cov(X_scaled.T)
    eig_values, eig_vectors = linalg.eig(cov)
    eig = zip(eig_values,eig_vectors)

    for k in eig:
        print(k)
