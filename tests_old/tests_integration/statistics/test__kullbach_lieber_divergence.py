import pytest


import numpy as np
import scipy.stats
from scipy.stats import norm
from scipy.stats import multivariate_normal
from scipy.stats import gaussian_kde
from pypospack.statistics import kullbach_lieber_divergence


# testing for normal distribution
def test__kld_calculation_1d_kde():
    n_samples_normal = 1000
    n_samples_kde = 1000
    rv_norm = norm(0,1)
    X_norm = rv_norm.rvs(size=1000)
    rv_kde_1 = gaussian_kde(X_norm)
    X_kde = rv_kde_1.resample(size=1000)
    rv_kde_2 = gaussian_kde(X_kde)
    kld = kullbach_lieber_divergence(rv_kde_1,rv_kde_2,1000)

    assert type(kld)==tuple
    assert kld[0]>0
    assert kld[0]>0

def dev__kld_calculation_1d_kde():
    n_samples_normal = 1000
    n_samples_kde = 1000
    rv_norm = norm(0,1)
    X_norm = rv_norm.rvs(size=1000)
    rv_kde_1 = gaussian_kde(X_norm)
    X_kde = rv_kde_1.resample(size=1000)
    rv_kde_2 = gaussian_kde(X_kde)
    kld = kullbach_lieber_divergence(rv_kde_1,rv_kde_2,1000)

    print(kld)

if __name__ == "__main__":
    dev__kld_calculation_1d_kde()
