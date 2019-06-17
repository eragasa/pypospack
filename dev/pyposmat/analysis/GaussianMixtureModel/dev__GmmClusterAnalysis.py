import os

import numpy as np
from numpy import linalg
import pandas as pd

import pypospack.utils
from pypospack.pyposmat.data import (PyposmatConfigurationFile, 
                                     PyposmatDataFile)
from gmm_analysis import GmmClusterAnalysis

# pypospack root directory
pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()

config_fn = os.path.join(
        pypospack_root_dir,
        'data','Si__sw__data','pareto_optimization_unconstrained',
        'pyposmat.config.in')
data_fn = os.path.join(
        pypospack_root_dir,
        'data','Si__sw__data','pareto_optimization_unconstrained',
        'pyposmat.kde.20.out')

n_components = 10
o = GmmClusterAnalysis(
        pyposmat_configuration=config_fn,
        pyposmat_data=data_fn,
        names='all',
        n_components=n_components)

# print(o.model.weights_)
# print(o.model.means_)
# print(o.model.covariances_)
# print(o.names)

assert len(set(o.data.df['cluster_id'].values)) == n_components

from scipy import stats
from sklearn.mixture import GaussianMixture 
qoi_names = o.configuration.qoi_names
qoi_indices = [o.names.index(k) for k in qoi_names]
qoi_target_values = np.array([o.configuration.qoi_targets[k] for k in qoi_names])

cluster_scores = []
for i in range(n_components):
    mean = o.model.means_[i][qoi_indices]
    covar = o.model.covariances_[i][np.ix_(qoi_indices,qoi_indices)]
    score = []
    for j in range(len(qoi_indices)):
        v = stats.norm.cdf(qoi_target_values[j],float(mean[j]),float(covar[j,j]))
        if v > .5:
            v = 1 - v
        score.append(v)
    ttl_score = 1
    for v in score:
        ttl_score = ttl_score * (1 - v)
    ttl_score = 1 - ttl_score
    score.append(ttl_score)
    cluster_scores.append(score)

if True:
    for score in cluster_scores:
        print(score)

#-----------------
# cluster info
#-------------
if True:
    
    header_row = ['cluster_id','weight','N']
    rows = []
    for i in range(n_components):
        row = [
                i,
                o.model.weights_[i],
                o.data.df.loc[o.data.df['cluster_id'] == i].shape[0]]
        rows.append(row)

if True:
    print(header_row)
    #str_out = ",".join(header_row) + "\n"
    for row in rows:
        print(row)
        #str_out += ",".join(row) + "\n"

#----
# parameter table
#----
if True:
    param_names =  o.configuration.free_parameter_names
    param_means = []
    param_stds = []
    for i in range(n_components):
        means = o.model.means_[i]
        covar = o.model.covariances_[i]
        param_means.append([
                o.model.means_[i][o.names.index(k)] 
                for k in param_names])
        param_stds.append([
                np.sqrt(o.model.covariances_[i][o.names.index(k),o.names.index(k)])
                for k in param_names])

if True:
    print(param_names)
    for i in range(n_components):
        print(param_means[i])
        print(param_stds[i])

#----------
# qoi table
#----------
if True:
    qoi_names =  o.configuration.qoi_names
    qoi_means = []
    qoi_stds = []
    for i in range(n_components):
        means = o.model.means_[i]
        covar = o.model.covariances_[i]
        qoi_means.append([
                o.model.means_[i][o.names.index(k)] 
                for k in qoi_names])
        qoi_stds.append([
                np.sqrt(o.model.covariances_[i][o.names.index(k),o.names.index(k)])
                for k in qoi_names])

if True:
    print(qoi_names)
    print(qoi_target_values)
    for i in range(n_components):
        print(qoi_means[i])
        print(qoi_stds[i])

#----------------
# best potentials 
#----------------
best_sim_ids = []
for i in range(n_components):
    qoi_names = o.configuration.qoi_names
    qoi_indices = [o.names.index(k) for k in qoi_names]
    mean = o.model.means_[i][qoi_indices]
    covar = o.model.covariances_[i][np.ix_(qoi_indices,qoi_indices)]
    df = o.data.df.loc[o.data.df['cluster_id']==i]
    X = o.data.df[qoi_names].loc[o.data.df['cluster_id']==i]
    df['score'] = stats.multivariate_normal.pdf(
            X,
            mean,
            covar)

    best_sim_ids += list(df.nlargest(2,'score')['sim_id'].values)

data = PyposmatDataFile()
data.read(filename=data_fn)
data.df['cluster_id'] = o.data.df['cluster_id']
data.df = data.df.loc[o.data.df['sim_id'].isin(best_sim_ids)]
data.write(filename='gmm_best_mle.out')

#---------------
#
#---------------
data = PyposmatDataFile()
data.read(filename=data_fn)
data.df['cluster_id'] = o.data.df['cluster_id']
data.df = data.df.loc[o.data.df['sim_id'].isin(best_sim_ids)]
data.write(filename='gmm_best_potentials.out')
qoi_names = o.configuration.qoi_names
qoi_targets = o.configuration.qoi_targets
for iqn,qn in enumerate(qoi_names):
    en = "{}.err".format(qn)
    nen = "{}.nerr".format(qn)
    q = qoi_targets[qn]
    data.df[nen] = data.df[qn]/q-1
normederr_names = ['{}.nerr'.format(q) for q in qoi_names]
data.df['score'] = np.sqrt(np.square(data.df[normederr_names]).sum(axis=1))

best_sim_ids = []
for i in range(n_components):
    best_sim_ids += list(data.df.loc[data.df['cluster_id'] == i].nsmallest(2,'score')['sim_id'].values)
data.df = data.df.loc[data.df['sim_id'].isin(best_sim_ids)]
data.write(filename='gmm_best_distance.out')

