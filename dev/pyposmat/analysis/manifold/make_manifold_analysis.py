import os
from collections import OrderedDict
import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np
import pandas as pd

from sklearn.mixture import GaussianMixture
from sklearn import preprocessing
from sklearn import manifold



if __name__ == "__main__":
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

    o_config = PyposmatConfigurationFile()
    o_config.read(filename=config_fn)
    
    o_data = PyposmatDataFile()
    o_data.read(filename=data_fn)

    manifold_learn_config = OrderedDict()
    manifold_learn_config['manifold_type'] = 'tsne'
    manifold_learn_config['pypospack_config_fn'] = config_fn 
    manifold_learn_config['pypospack_data_fn'] =  data_fn


    fig,ax = plt.subplots(1,3)


    if manifold_learn_config['manifold_type'] == 'mds':
        manifold['config'] = OrderedDict()
        manifold['config']['n_components'] = 2
        manifold['config']['max_iter'] = 1000
        manifold['config']['n_init'] = 1

        print('parameter_analysis')
        X_param = o_data.df[o_config.free_parameter_names]
        X_param_scaled = preprocessing.scale(X_param)
        mds_n_components = 2
        mds_param = MDS(n_components = 2,max_iter=1000, n_init=1)
        Y_param = mds_param.fit_transform(X_param_scaled)
        ax[0].scatter(Y_param[:,0],Y_param[:,1],s=1,c='black')

        print('qoi_analysis')
        X_qoi = o_data.df[o_config.qoi_names]
        X_qoi_scaled = preprocessing.scale(X_param)
        mds_n_components = 2
        mds_qoi = MDS(n_components = 2,max_iter=1000, n_init=1)
        Y_qoi = mds_qoi.fit_transform(X_qoi_scaled)
        ax[1].scatter(Y_qoi[:,0],Y_qoi[:,1],s=1,c='black')

        print('parameter+qoi')
        X_all = o_data.df[o_config.free_parameter_names + o_config.qoi_names]
        X_all_scaled = preprocessing.scale(X_param)
        mds_n_components = 2
        mds_all = MDS(n_components = 2,max_iter=1000, n_init=1)
        Y_all = mds_all.fit_transform(X_all_scaled)
        ax[2].scatter(Y_all[:,0],Y_all[:,1],s=1,c='black')

    elif manifold_learn_config['manifold_type'] == 'tsne':
        manifold_learn_config['config'] = OrderedDict()
        manifold_learn_config['config']['n_components'] = 2
        manifold_learn_config['config']['init'] = 'pca'
        manifold_learn_config['config']['random_state'] = 0

        print('parameter_analysis')
        X_param = o_data.df[o_config.free_parameter_names]
        Y_param = manifold.TSNE(**manifold_learn_config['config']).fit_transform(X_param)
        print(X_param.shape)
        print(Y_param.shape)
        ax[0].scatter(Y_param[:,0],Y_param[:,1],s=1,c='black')

        print('qoi_analysis')
        X_qoi = o_data.df[o_config.qoi_names]
        Y_qoi = manifold.TSNE(**manifold_learn_config['config']).fit_transform(X_qoi)
        print(X_qoi.shape)
        print(Y_qoi.shape)
        ax[1].scatter(Y_qoi[:,0],Y_qoi[:,1],s=1,c='black')

        print('parameter+qoi')
        X_all = o_data.df[o_config.free_parameter_names + o_config.qoi_names]
        Y_all = manifold.TSNE(**manifold_learn_config['config']).fit_transform(X_all)
        print(X_all.shape)
        print(Y_all.shape)
        ax[2].scatter(Y_all[:,0],Y_all[:,1],s=1,c='black')

    for i in range(len(ax)):
        ax[i].set_aspect('equal',adjustable='box')
        ax[i].xaxis.set_major_formatter(NullFormatter())
        ax[i].yaxis.set_major_formatter(NullFormatter())
    fig.tight_layout()
    plt.show()



