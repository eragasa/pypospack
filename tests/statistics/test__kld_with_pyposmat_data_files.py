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
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile

test_case = OrderedDict()
test_case['Si__sw'] = OrderedDict()
test_case['Si__sw']['data_directory'] = os.path.join(
        pypospack.utils.get_pypospack_root_directory(),
        'examples/Si__sw/data'
)
test_case['Si__sw']['configuration_fn'] = os.path.join(test_case['Si__sw']['data_directory'],'pyposmat.config.in')
test_case['Si__sw']['kde_files'] = [f for f in os.listdir(test_case['Si__sw']['data_directory']) if f.startswith('pyposmat.kde')]
test_case['Si__sw']['results_files'] = [f for f in os.listdir(test_case['Si__sw']['data_directory']) if f.startswith('pyposmat.results')]

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
    print(test_case['Si__sw']['data_directory'])
    print(test_case['Si__sw']['results_files'].sort())
    print(test_case['Si__sw']['kde_files'].sort())

    n_results_files = len(test_case['Si__sw']['results_files'])
    n_kde_files = len(test_case['Si__sw']['kde_files'])

    configuration_fn = test_case['Si__sw']['configuration_fn']
    configuration = PyposmatConfigurationFile()
    configuration.read(filename=configuration_fn)
    free_parameter_names = configuration.free_parameter_names
    print('free_parameter_names:{}'.format(free_parameter_names))
    
    for i in range(n_kde_files):
        if i > 0:
            print(80*'-')
            print('i_iteration:{}'.format(i))

            data_directory = test_case['Si__sw']['data_directory']
            kde_file_fn_1 = os.path.join(data_directory,'pyposmat.kde.{}.out'.format(i))
            kde_file_fn_2 = os.path.join(data_directory,'pyposmat.kde.{}.out'.format(i+1))

            print('kde_file_fn_1:{}'.format(kde_file_fn_1))
            print('kde_file_fn_2:{}'.format(kde_file_fn_2))
            

            kde_file_1 = PyposmatDataFile()
            kde_file_1.read(filename=kde_file_fn_1)
            kde_rv_1 = gaussian_kde(kde_file_1.df[free_parameter_names].T)

            kde_file_2 = PyposmatDataFile()
            kde_file_2.read(filename=kde_file_fn_2)
            kde_rv_2 = gaussian_kde(kde_file_2.df[free_parameter_names].T)

            kld = kullbach_lieber_divergence(kde_rv_1,kde_rv_2,1000)
            print(kld)
    dev__kld_calculation_1d_kde()
