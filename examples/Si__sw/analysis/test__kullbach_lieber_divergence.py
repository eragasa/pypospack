import pytest

import os
from collections import OrderedDict
import numpy as np
import scipy.stats
from scipy.stats import norm
from scipy.stats import multivariate_normal
from scipy.stats import gaussian_kde

import pypospack.utils
from pypospack.statistics import kullbach_lieber_divergence
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

def calculate_kld(config, data_1_fn,data_2_fn,n_samples=10000):
    assert isinstance(config,PyposmatConfigurationFile) or isinstance(config,str)
    assert isinstance(data_1_fn,str)
    assert isinstance(data_2_fn,str)
    assert isinstance(n_samples,int)

    assert os.path.isfile(data_1_fn)
    assert os.path.isfile(data_1_fn)

    if isinstance(config,PyposmatConfigurationFile):
        o_config = config
    elif isinstance(config, str):
        o_config = PyposmatConfigurationFile()
        o_config.read(filename=config)

    # initialize
    kld = OrderedDict()
    
    data_1 = PyposmatDataFile()
    data_1.read(filename=data_1_fn)

    data_2 = PyposmatDataFile()
    data_2.read(filename=data_2_fn)

    kde_1 = gaussian_kde(data_1.df[o_config.free_parameter_names].T)
    kde_2 = gaussian_kde(data_2.df[o_config.free_parameter_names].T)
    kld['param'] = kullbach_lieber_divergence(kde_1,kde_2,n_samples)[0]

    kde_1 = gaussian_kde(data_1.df[o_config.qoi_names].T)
    kde_2 = gaussian_kde(data_2.df[o_config.qoi_names].T)
    kld['qoi'] = kullbach_lieber_divergence(kde_1,kde_2,n_samples)[0]
    
    return kld


data_directory = os.path.join(
        pypospack.utils.get_pypospack_root_directory(),
        'data','Si__sw__data','pareto_optimization_unconstrained')

config_fn = os.path.join(data_directory,'pyposmat.config.in')

o_config = PyposmatConfigurationFile()
o_config.read(filename=config_fn)

for i in range(o_config.n_iterations):
    if i == 0:
        kld_results = None
        kld_kde = None
        kld_filter = calculate_kld(
                config=o_config,
                data_1_fn=os.path.join(data_directory,'pyposmat.results.{}.out'.format(i)),
                data_2_fn=os.path.join(data_directory,'pyposmat.kde.{}.out'.format(i+1)),
                n_samples=10000)
    else:
        kld_results = calculate_kld(
                config=o_config,
                data_1_fn=os.path.join(data_directory,'pyposmat.results.{}.out'.format(i-1)),
                data_2_fn=os.path.join(data_directory,'pyposmat.results.{}.out'.format(i)),
                n_samples=10000)
        kld_kde = calculate_kld(
                config=o_config,
                data_1_fn=os.path.join(data_directory,'pyposmat.kde.{}.out'.format(i)),
                data_2_fn=os.path.join(data_directory,'pyposmat.kde.{}.out'.format(i+1)),
                n_samples=10000)
        kld_filter = calculate_kld(
                config=o_config,
                data_1_fn=os.path.join(data_directory,'pyposmat.results.{}.out'.format(i)),
                data_2_fn=os.path.join(data_directory,'pyposmat.kde.{}.out'.format(i+1)),
                n_samples=10000)
    print("{} {} {} {}".format(i,kld_results,kld_kde,kld_results,kld_filter))

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
    pass
    # dev__kld_calculation_1d_kde()
