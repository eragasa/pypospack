import os
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np
import pandas as pd

from sklearn.mixture import GaussianMixture
from sklearn import preprocessing
from sklearn import manifold

class Manifold(object):
    
    def __init__(self,pyposmat_configuration,pyposmat_data,manifold_config=None):
        self.configuration = None
        self.data = None
        self.manifold_configuration = None
        self.initialize_configuration(pyposmat_configuration=pyposmat_configuration)
        self.initialize_data(pyposmat_data=pyposmat_data)
        self.initialize_manifold_config(manifold_config=manifold_config)

    def initialize_configuration(self,pyposmat_configuration):
        if isinstance(pyposmat_configuration,PyposmatConfigurationFile):
            self.configuration = pyposmat_configuration
        elif isinstance(pyposmat_configuration,str):
            self.configuration = PyposmatConfigurationFile()
            self.configuration.read(filename=pyposmat_configuration)
        else:
            raise TypeError('pyposmat_configuration must be either a path or PyposmatConfigurationFile')

    def initialize_data(self,pyposmat_data):
        if isinstance(pyposmat_data, PyposmatDataFile):
            self.data = pyposmat_data
        elif isinstance(pyposmat_data, str):
            self.data = PyposmatDataFile
            self.data.read(filename=pyposmat_data)
        else:
            raise TypeError('pyposmat_data must either be a path or a PyposmatDataFile')

    def initialize_manifold_configuration(self,manifold_configuration=None):
        if manifold_configuration is None:
            self.manifold_configuration = None
        else:
            raise NotImplementedError

    def learn_manifold(self,names,scaling_type):
        raise NotImplementedError

    def scale_data(self,X,scaling_type):
        if scaling_type='standard':
            X_scaled = preprocessing.scale(X)
        elif scaling_type='none':
            X_scaled = X
        else:
            raise ValueError('unknown scaling type')

class ManifoldLearningFactory(object):
    pass

class ManifoldPlotter(object):
    pass

class MdsManifold(ManifoldAnalysis):

    manifold_type = 'MDS'

    def __init__(self,pyposmat_configuration,pyposmat_data,manifold_config=None):
        Manifold.__init__(self,
                pyposmat_configuration=pyposmat_configuration,
                pyposmat_data=pyposmat_data,
                manifold_config=manifold_config
                )

    def initialize_manifold_configuration(self,manifold_configuration=None):
        DEFAULT_MDS_MANIFOLD_CONFIGURATION = OrderedDict()
        DEFAULT_MDS_MANIFOLD_CONFIGURATION['n_components'] = 2
        DEFAULT_MDS_MANIFOLD_CONFIGURATION['init'] = 'pca'
        DEFAULT_MDS_MANIFOLD_CONFIGURATION['random_state'] = 0

        if manifold_configuration is None:
            self.manifold_configuration = DEFAULT_TSNE_MANIFOLD_CONFIGURATION
        elif isinstance(manifold_configuration,dict):
            self.manifold_configuration = manifold_configuration

            for k,v in DEFAULT_MDS_MANIFOLD_CONFIGURATION.items():
                if k not in self.manifold_configuration:
                    self.manifold_configuration[k] = v

    def learn_manifold(self,names,scaling_type='standard'):
        if isinstance(names,list):
            names_ = names
        elif names == 'qois':
            names_ = self.configuration.qoi_names
        elif names == 'free_parameters':
            names_ = self.configuration.free_parmaeter_names
        elif names == 'all':
            names_ = self.configuration.qoi_names \
                    + self.configuration.free_parameter_names

        X = self.data.df[names_]
        X_scaled = get_scaled_data(X,scaling_type=scaling_type)
        self.manifold = manifold.MDS(**self.manifold_configuration)
        Y = mds.fit_transform(X_scaled)
        
        for i in self.manifold_configuration['n_components']:
            self.data.df['MDS_{}'.format(str(i))] = Y[:,i]

    def MDS(self,i):
        return self.data.df['MDS_{}'.format(str(i))]

    def scale_data(self,X,scaling_type):
        if scaling_type='standard':
            X_scaled = preprocessing.scale(X)
        elif scaling_type='none':
            X_scaled = X
        else:
            raise ValueError('unknown scaling type')

class TsneManifold(ManifoldAnaylsis):
   
    manifold_type = 'TSNE'
    def __init__(self,pyposmat_configuration,pyposmat_data,manifold_config=None):
        Manifold.__init__(self,
                pyposmat_configuration=pyposmat_configuration,
                pyposmat_data=pyposmat_data,
                manifold_config=manifold_config
                )

    def initialize_manifold_configuration(self,manifold_configuration=None):
        DEFAULT_TSNE_MANIFOLD_CONFIGURATION = OrderedDict()
        DEFAULT_TSNE_MANIFOLD_CONFIGURATION['n_components'] = 2
        DEFAULT_TSNE_MANIFOLD_CONFIGURATION['init'] = 'pca'
        DEFAULT_TSNE_MANIFOLD_CONFIGURATION['random_state'] = 0
        DEFAULT_TSNE_MANIFOLD_CONFIGURATION['perplexity'] = 30
        DEFAULT_TSNE_MANIFOLD_CONFIGURATION['learning_rate'] = 200
        DEFAULT_TSNE_MANIFOLD_CONFIGURATION['n_iter'] = 1000
        DEFAULT_TSNE_MANIFOLD_CONFIGURATION['n_iter_without_progress'] = 300
        DEFAULT_TSNE_MANIFOLD_CONFIGURATION['min_grad_norm'] = 1e-7
        DEFAULT_TSNE_MANIFOLD_CONFIGURATION['metric'] = 'euclidean'
        DEFAULT_TSNE_MANIFOLD_CONFIGURATION['init'] = 'pca'

        if manifold_configuration is None:
            self.manifold_configuration = DEFAULT_TSNE_MANIFOLD_CONFIGURATION
        elif isinstance(manifold_configuration,dict):
            self.manifold_configuration = manifold_configuration

            for k,v in DEFAULT_TSNE_MANIFOLD_CONFIGURATION.items():
                if k not in self.manifold_configuration:
                    self.manifold_configuration[k] = v

    def get_manifold_properties(self):
        d = OrdereDict()
        d['kl_divergence'] = self.manifold.kl_divergence
        d['n_iter'] = self.manifold.n_iter_
        
        return d
    def learn_manifold(self,names,scaling_type=None):
        if isinstance(names,list):
            names_ = names
        elif names == 'qois':
            names_ = self.configuration.qoi_names
        elif names == 'free_parameters':
            names_ = self.configuration.free_parmaeter_names
        elif names == 'all':
            names_ = self.configuration.qoi_names \
                    + self.configuration.free_parameter_names

        if scaling_type is None:
            pass
        elif scaling_type == 'pca':
            self.manifold_configuration['init'] == 'pca'

        X = self.data.df[names_]
        self.manifold = manifold.TSNE(**self.manifold_configuration)
        Y = mds.fit_transform(X_scaled)
        
        for i in self.manifold_configuration['n_components']:
            data_name = 'TSNE_{}'.format(str(i))
            data_type = 'TSNE'
            self.data.df[data_name] = Y[:,i]

    def TSNE(self,i):
        return self.data.df['TSNE_{}'.format(str(i))]

if __name__ == "__main__":
    import pypospack.utils
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



