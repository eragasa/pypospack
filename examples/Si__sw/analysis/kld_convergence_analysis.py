import pytest

import os
from collections import OrderedDict
import numpy as np
from numpy import linalg
import scipy.stats
from scipy.stats import norm
from scipy.stats import multivariate_normal
from scipy.stats import gaussian_kde

import pypospack.utils
from pypospack.statistics import kullbach_lieber_divergence
from pypospack.statistics import GaussianKde
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile


def covariance_analysis(data_fn,names):
    assert isinstance(data_fn,str)
    assert isinstance(names,list)

    data = PyposmatDataFile()
    data.read(filename=data_fn)

    cov_matrix = np.cov(data.df[names].T)
    w,v = linalg.eig(cov_matrix)
    print("eigenvalues:\n",w)
    print("eigenvectors:\n",v)

def calculate_kld(data_1_fn,data_2_fn,names,n_samples=2000):
    assert isinstance(data_1_fn,str)
    assert isinstance(data_2_fn,str)
    assert isinstance(n_samples,int)

    assert os.path.isfile(data_1_fn)
    assert os.path.isfile(data_1_fn)

    data_1 = PyposmatDataFile()
    data_1.read(filename=data_1_fn)

    data_2 = PyposmatDataFile()
    data_2.read(filename=data_2_fn)

    w1,v1 = linalg.eig(np.cov(data_1.df[names].T))
    w2,v2 = linalg.eig(np.cov(data_2.df[names].T))
  
    cov1_ill_conditioned = any([k < 0 for k in w1.tolist()])
    cov2_ill_conditioned = any([k < 0 for k in w2.tolist()])

    any_ill_conditioned = any([cov1_ill_conditioned,cov2_ill_conditioned])

    if any_ill_conditioned:
        kde_1 = GaussianKde(data_1.df[names].T)
        kde_2 = GaussianKde(data_2.df[names].T)
    else:
        kde_1 = gaussian_kde(data_1.df[names].T)
        kde_2 = gaussian_kde(data_2.df[names].T)
    
    kld = kullbach_lieber_divergence(kde_1,kde_2,n_samples)
    return kld

def calculate_kld_parameters(config,data_directory,kld_param_fn='pyposmat.kld_param.out'):
    assert isinstance(config,str) or isinstance(config,PyposmatConfigurationFile)
    assert os.path.isdir(data_directory)
    assert isinstance(kld_param_fn,str)

    # process the the configuration argument, the configuration argument has two
    # options for processing
    # (1) if config is a str, the config is assumed to be a path to the
    #     the path to the configuration file, and o_config is initialized from it
    # (2) if config is PyposmatConfigurationFile object, then o_config is set to it
    if isinstance(config,str):
        o_config = PyposmatConfigurationFile()
        o_config.read(filename=config)
    else:
        assert isinstance(config,PyposmatConfigurationFile)
        o_config = config

    kld = OrderedDict()
    for i in range(o_config.n_iterations):
        kld[i] = OrderedDict()
        if i == 0:
            kld[i]['results'] = None
            kld[i]['kde'] = None

            kld[i]['filter'] = calculate_kld(
                data_1_fn=os.path.join(data_directory,'pyposmat.results.{}.out'.format(i)),
                data_2_fn=os.path.join(data_directory,'pyposmat.kde.{}.out'.format(i+1)),
                names = o_config.free_parameter_names,
                n_samples=n_samples)
        else:
            kld[i]['results'] = calculate_kld(
                data_1_fn=os.path.join(data_directory,'pyposmat.results.{}.out'.format(i-1)),
                data_2_fn=os.path.join(data_directory,'pyposmat.results.{}.out'.format(i)),
                names = o_config.free_parameter_names,
                n_samples=n_samples)
            kld[i]['kde'] = calculate_kld(
                data_1_fn=os.path.join(data_directory,'pyposmat.kde.{}.out'.format(i)),
                data_2_fn=os.path.join(data_directory,'pyposmat.kde.{}.out'.format(i+1)),
                names = o_config.free_parameter_names,
                n_samples=n_samples)
            kld[i]['filter'] = calculate_kld(
                data_1_fn=os.path.join(data_directory,'pyposmat.results.{}.out'.format(i)),
                data_2_fn=os.path.join(data_directory,'pyposmat.kde.{}.out'.format(i+1)),
                names = o_config.free_parameter_names,
                n_samples=n_samples)
        print(i)

    # write out the kld_parameters file
    with open(kld_param_fn,'w') as f:
        f.write(",".join(['iteration','results','kde','filter'])+"\n")

        for kld_iteration,kld_row in kld.items():
            s_list = []
            s_list.append(kld_iteration)
            for k in ['results','kde','filter']:
                if kld_row[k] is None:
                    s_list.append(float('NaN'))
                else:
                    s_list.append(kld_row[k][0])
            f.write(",".join([str(s) for s in s_list])+"\n")

def calculate_kld_qois(config,data_directory,kld_qoi_fn='pyposmat.kld_qoi.out'):
    assert isinstance(config,str) or isinstance(config,PyposmatConfigurationFile)
    assert os.path.isdir(data_directory)
    assert isinstance(kld_param_fn,str)
  
    if isinstance(config,str):
        o_config = PyposmatConfigurationFile()
        o_config.read(filename=config)
    else:
        assert isinstance(config,PyposmatConfigurationFile)
        o_config = config

    kld = OrderedDict()
    s = [",".join(['iteration','results','kde','filter'])]
    print(s)
    for i in range(o_config.n_iterations):
        kld[i] = OrderedDict()
        if i == 0:
            kld[i]['results'] = None
            kld[i]['kde'] = None
            kld[i]['filter'] = calculate_kld(
                data_1_fn=os.path.join(data_directory,'pyposmat.results.{}.out'.format(i)),
                data_2_fn=os.path.join(data_directory,'pyposmat.kde.{}.out'.format(i+1)),
                names = [k for k in o_config.qoi_names if k != "Si_dia.B"],
                n_samples=n_samples)
        else:
            data_1_path = os.path.join(data_directory,'pyposmat.results.{}.out'.format(i-1))
            data_2_path = os.path.join(data_directory,'pyposmat.results.{}.out'.format(i))
            
            o_data_1 = PyposmatDataFile()
            o_data_1.read(filename=data_1_path)

            kld[i]['results'] = calculate_kld(
                data_1_fn=os.path.join(data_directory,'pyposmat.results.{}.out'.format(i-1)),
                data_2_fn=os.path.join(data_directory,'pyposmat.results.{}.out'.format(i)),
                names = o_config.qoi_names,
                n_samples=n_samples)
            kld[i]['kde'] = calculate_kld(
                data_1_fn=os.path.join(data_directory,'pyposmat.kde.{}.out'.format(i)),
                data_2_fn=os.path.join(data_directory,'pyposmat.kde.{}.out'.format(i+1)),
                names = o_config.qoi_names,
                n_samples=n_samples)
            kld[i]['filter'] = calculate_kld(
                data_1_fn=os.path.join(data_directory,'pyposmat.results.{}.out'.format(i)),
                data_2_fn=os.path.join(data_directory,'pyposmat.kde.{}.out'.format(i+1)),
                names = o_config.qoi_names,
                n_samples=n_samples)
        
        s_row = [i]
        for k in ['results','kde','filter']:
            if kld[i][k] is None:
                s_row += [None]
            else:
                s_row+= [kld[i][k][0]]
        s = ",".join([str(v) for v in s_row])
        print(s)

    if not os.path.isfile(kld_qoi_fn):
        with open(kld_param_fn,'w') as f:
            f.write(",".join(['iteration','results','kde','filter'])+"\n")

            for kld_iteration,kld_row in kld.items():
                s_list = []
                s_list.append(kld_iteration)
                for k in ['results','kde','filter']:
                    if kld_row[k] is None:
                        s_list.append(float('NaN'))
                    else:
                        s_list.append(kld_row[k][0])
                f.write(",".join([str(s) for s in s_list])+"\n")

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
    print(80*'-')
    print('{:^80}'.format('dev__kld_calculation_1d_kde'))
    print(80*'-')

    n_samples_normal = 1000
    n_samples_kde = 1000
    rv_norm = norm(0,1)
    X_norm = rv_norm.rvs(size=1000)
    rv_kde_1 = gaussian_kde(X_norm)
    X_kde = rv_kde_1.resample(size=1000)
    rv_kde_2 = gaussian_kde(X_kde)

    assert isinstance(rv_kde_1,gaussian_kde)
    assert isinstance(rv_kde_1,gaussian_kde)
    kld = kullbach_lieber_divergence(rv_kde_1,rv_kde_2,1000)

    print(kld)
    return kld

if __name__ == "__main__":
    kld_output_path = "kld_analysis"
    try:
        os.mkdir(kld_output_path)
    except FileExistsError as e:
        if os.path.isdir(kld_output_path):
            pass
        else:
            raise
    except Exception as e:
        raise

    dev__kld_calculation_1d_kde()
    n_samples = 10000
    data_directory = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'data','Si__sw__data','pareto_optimization_unconstrained')

    config_fn = os.path.join(data_directory,'pyposmat.config.in')
    o_config = PyposmatConfigurationFile()
    o_config.read(filename=config_fn)

    kld_param_fn = 'pyposmat.kld_param.out'
    kld_qoi_fn = 'pyposmat.kld_qoi.out'

    if not os.path.isfile(kld_param_fn):
        kld_param_fn = os.path.join(kld_output_path,'pyposmat.kld_param.out')
        calculate_kld_parameters(
                config=o_config,
                data_directory=data_directory,
                kld_param_fn=kld_param_fn)

    if not os.path.isfile(kld_qoi_fn):
        kld_qoi_fn = os.path.join(kld_output_path,'pyposmat.kld_qoi.out')
        calculate_kld_qois(
                config=o_config,
                data_directory=data_directory,
                kld_qoi_fn=kld_qoi_fn)
    
