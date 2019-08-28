import pytest

import numpy as np
from scipy import stats

from pypospack.statistics import GaussianKde

def measure(n):
    "Measurement model, return two coupled measurements."
    m1 = np.random.normal(size=n)
    m2 = np.random.normal(scale=0.5, size=n)
    return m1+m2, m1-m2

@pytest.fixture()
def resource_measure():
   m1,m2 = measure(2000)
   X = np.vstack([m1,m2])
   return X

def test____init__(resource_measure):
   X = resource_measure
   n, d = X.shape

   o = GaussianKde(X)
   assert o.n == n
   assert o.d == d

def test____init____scaler_means(resource_measure):
   X = resource_measure
   n, d = X.shape

   o = GaussianKde(X)
   assert o.n == n
   assert o.d == d

   for i in range(d):
       assert X[i,:].mean() - o.scaler.mean_[i] < 1e-6

def test____init____scaler_std(resource_measure):
   X = resource_measure
   n, d = X.shape

   o = GaussianKde(X)
   assert o.n == n
   assert o.d == d

   for i in range(d):
       assert X[i,:].std() - np.sqrt(o.scaler.var_[i]) < 1e-6

@pytest.mark.parametrize('bw_method',[(k) for k in GaussianKde.bw_method_options])
def test____init____bw_methods(resource_measure,bw_method):
    X = resource_measure
    n, d = X.shape

    o = GaussianKde(X,bw_method=bw_method)
    assert o.n == n
    assert o.d == d

def test__evaluate(resource_measure):
    X = resource_measure
    n, d = X.shape

    o = GaussianKde(X)

    
if __name__ == "__main__":

    m1,m2 = measure(2000)
    X = np.vstack([m1,m2])
    stats.gaussian_kde(X,bw_method='scott', weights=[1,2,3])

