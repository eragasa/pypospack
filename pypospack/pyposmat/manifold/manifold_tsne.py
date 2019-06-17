# -*- coding: utf-8 -*-
"""Implementation of TSNE manifold class

This module provides the abstract class implementation for manifold learning
this wrapper is primarily a wrapper from the scikit-learn module.
"""
from collections import OrderedDict
from pypospack.pyposmat.manifold import Manifold

from sklearn import manifold

class TsneManifold(Manifold):
   
    manifold_type = 'TSNE'

    DEFAULT_MANIFOLD_CONFIGURATION = OrderedDict(
            [
                ('n_components', 2),
                ('init',' pca'),
                ('random_state', 0),
                ('perplexity', 50),
                ('learning_rate', 10),
                ('n_iter', 5000),
                ('n_iter_without_progress', 500),
                ('min_grad_norm', 1e-8),
                ('metric', 'euclidean'),
                ('init', 'pca'),
                ('verbose', 2)
            ]
        )

    def __init__(self,
                 pyposmat_configuration,
                 pyposmat_data,
                 manifold_configuration=None):

        Manifold.__init__(self,
                pyposmat_configuration=pyposmat_configuration,
                pyposmat_data=pyposmat_data,
                manifold_configuration=manifold_configuration
                )

    def initialize_manifold_configuration(self, configuration=None):
        """

        Overrides the default behavior of the abstract Manifold class

        """
        DEFAULT_TSNE_MANIFOLD_CONFIGURATION = OrderedDict()
        DEFAULT_TSNE_MANIFOLD_CONFIGURATION['n_components'] = 2
        DEFAULT_TSNE_MANIFOLD_CONFIGURATION['init'] = 'pca'
        DEFAULT_TSNE_MANIFOLD_CONFIGURATION['random_state'] = 0
        DEFAULT_TSNE_MANIFOLD_CONFIGURATION['perplexity'] = 50
        DEFAULT_TSNE_MANIFOLD_CONFIGURATION['learning_rate'] = 10
        DEFAULT_TSNE_MANIFOLD_CONFIGURATION['n_iter'] = 5000
        DEFAULT_TSNE_MANIFOLD_CONFIGURATION['n_iter_without_progress'] = 500
        DEFAULT_TSNE_MANIFOLD_CONFIGURATION['min_grad_norm'] = 1e-8
        DEFAULT_TSNE_MANIFOLD_CONFIGURATION['metric'] = 'euclidean'
        DEFAULT_TSNE_MANIFOLD_CONFIGURATION['init'] = 'pca'
        DEFAULT_TSNE_MANIFOLD_CONFIGURATION['verbose'] = 2

        if configuration is None:
            self.manifold_configuration = DEFAULT_TSNE_MANIFOLD_CONFIGURATION
        elif isinstance(configuration,dict):
            self.manifold_configuration = configuration

            for k,v in DEFAULT_TSNE_MANIFOLD_CONFIGURATION.items():
                if k not in self.manifold_configuration:
                    self.manifold_configuration[k] = v

    def get_manifold_properties(self):
        d = OrderedDict()
        d['kl_divergence'] = self.manifold.kl_divergence_
        #d['n_iter'] = self.manifold.n_iter_
        
        return d

    def learn_manifold(self,names,scaling_type=None):
        if isinstance(names,list):
            names_ = names
        elif names == 'qois':
            names_ = self.configuration.qoi_names
        elif names == 'free_parameters':
            names_ = self.configuration.free_parameter_names
        elif names == 'all':
            names_ = self.configuration.qoi_names \
                    + self.configuration.free_parameter_names

        if scaling_type is None:
            pass
        elif scaling_type == 'pca':
            self.manifold_configuration['init'] == 'pca'

        X = self.data.df[names_]
        X_scaled = X
        self.manifold = manifold.TSNE(**self.manifold_configuration)
        Y = self.manifold.fit_transform(X_scaled)
        
        for i in range(self.manifold_configuration['n_components']):
            data_name = 'TSNE_{}'.format(str(i+1))
            data_type = 'TSNE'
            self.data.df[data_name] = Y[:,i]

    def TSNE(self,i):
        return self.data.df['TSNE_{}'.format(str(i))]

if __name__ == "__main__":
    pass
