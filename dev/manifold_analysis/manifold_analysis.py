import os
from copy import deepcopy
from collections import OrderedDict

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN
from sklearn.manifold import TSNE
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

class ManifoldAnalysis(object):

    def __init__(self,
                 configuration,
                 data,
                 n_components=2,
                 names=None,
                 output_path='pca_analysis'):
        self._initialize_configuration(configuration=configuration)
        self._initialize_data(data=data)
        self._initialize_names(names=names)
        self.output_path = output_path

        assert isinstance(n_components, int)
        self.n_components = n_components

        self.scaler = None
        self.pca = None
        self.cluster_ids = None

    def _initialize_configuration(self, configuration):
        if isinstance(configuration, str):
            assert os.path.isfile(configuration)
            self.configuration = PyposmatConfigurationFile()
            self.configuration.read(filename=configuration)
        elif isinstance(configuration, PyposmatConfigurationFile):
            self.configuration = configuration
        else:
            raise TypeError('configuration cannot be type:{}'.format(str(type(configuration))))

    def _initialize_data(self, data):
        if isinstance(data, str):
            assert os.path.isfile(data)
            self.data = PyposmatDataFile()
            self.data.read(filename=data)
        elif isinstance(data, PyposmatDataFile):
            self.data = deepcopy(data)
        else:
            raise TypeError('data cannot be type:{}'.format(str(type(data))))

        self.data.create_normalized_errors(
                normalize_type='by_qoi_target',
                qoi_targets=self.configuration.qoi_targets)

    def _initialize_names(self, names):
        if isinstance(names, list):
            self.names = list(names)
        elif isinstance(names, str):
            if names == 'qois':
                self.names = self.configuration.qoi_names
            elif names == 'parameters':
                self.names = self.configuration.parameter_names
            elif names == 'all':
                self.names = self.configuration.parameter_names + self.configuration.qoi_names
            else:
                raise TypeError
        else:
            raise TypeError

    def make_manifold(self, manifold_type='tsne', n_components=2, names='all'):
        if names == 'parameters':
            names_ = self.configuration.parameter_names
        elif names == 'qois':
            names_ = self.configuration.qoi_names
        elif names == 'all':
            names_ = self.configuration.parameter_names + self.configuration.qoi_names
        else:
            raise ValueError

        data_ = self.data.df[names_]

        if manifold_type == 'tsne':
            manifold = TSNE(n_components=n_components)
        
        manifold_components = manifold.fit_transform(data_)

        if names == 'parameters':
            for i in range(n_components):
                column_name = 'tsne_param_{}'.format(i+1)
                self.data.df[column_name] = manifold_components[:,i]
        elif names == 'qois':
            for i in range(n_components):
                self.data.df['tsne_qoi_{}'.format(i+1)] = manifold_components[:,i]
        elif names == 'all':
            for i in range(n_components):
                self.data.df['tsne_{}'.format(i+1)] = manifold_components[:,i]
        else:
            raise ValueError

    def plot_manifold_analysis(self):

        fig, ax = plt.subplots(1,3)

        print('learning parameter manifold...')
        self.make_manifold(names='parameters')
        print('learning qoi manifold...')
        self.make_manifold(names='qois')
        print('learning total manifold...')
        self.make_manifold(names='all')

        ax[0].scatter(self.data.df['tsne_param_1'],
                      self.data.df['tsne_param_2'],
                      s=1.)
        ax[0].set_xlabel('tsne_param_1')
        ax[0].set_ylabel('tsne_param_2')

        ax[1].scatter(self.data.df['tsne_qoi_1'],
                      self.data.df['tsne_qoi_2'],
                      s=1.)
        ax[1].set_xlabel('tsne_qoi_1')
        ax[1].set_ylabel('tsne_qoi_2')
        
        ax[2].scatter(self.data.df['tsne_1'],
                      self.data.df['tsne_2'],
                      s=1.)
        ax[2].set_xlabel('tsne_1')
        ax[2].set_ylabel('tsne_2')

        fig.tight_layout()
        plt.show()

    def plot_cluster_analysis(self,
                              cluster_type='dbscan',
                              names='all'):
        fig, ax = plt.subplots(1,3)

        print('learning parameter manifold...')
        self.make_manifold(names='parameters')
        print('learning qoi manifold...')
        self.make_manifold(names='qois')
        print('learning total manifold...')
        self.make_manifold(names='all')

        cluster_type = 'dbscan'
        if cluster_type == 'dbscan':
            dbscan_args = {
                    'eps':0.5
            }
            cluster = DBSCAN(**dbscan_args)
        elif cluster_type == 'optics':
            cluster = OPTICS()

        if names == 'parameters':
            cols = ['tsne_param_1', 'tsne_param_2']
            data_ =self.data.df[cols]
            self.data.df['cluster_id'] = cluster.fit_predict(data_)
        elif names == 'qois':
            cols = ['tsne_qoi_1', 'tsne_qoi_2']
            data_ =self.data.df[cols]
            self.data.df['cluster_id'] = cluster.fit_predict(data_)
        elif names == 'all':
            cols = ['tsne_1', 'tsne_2']
            data_ =self.data.df[cols]
            self.data.df['cluster_id'] = cluster.fit_predict(data_)
        else:
            raise ValueError

        self.cluster_ids = list(set(self.data.df['cluster_id'].values))
        self.n_clusters = len(self.cluster_ids)
        print(self.cluster_ids)
        for i in self.cluster_ids:
            ax[0].scatter(self.data.df['tsne_param_1'].loc[self.data.df['cluster_id'] == i],
                          self.data.df['tsne_param_2'].loc[self.data.df['cluster_id'] == i],
                          s=1.)
            ax[1].scatter(self.data.df['tsne_qoi_1'].loc[self.data.df['cluster_id'] == i],
                          self.data.df['tsne_qoi_2'].loc[self.data.df['cluster_id'] == i],
                          s=1.)
            ax[2].scatter(self.data.df['tsne_1'].loc[self.data.df['cluster_id'] == i],
                          self.data.df['tsne_2'].loc[self.data.df['cluster_id'] == i],
                          s=1.)

        ax[0].set_xlabel('tsne_param_1')
        ax[0].set_ylabel('tsne_param_2')
        ax[1].set_xlabel('tsne_qoi_1')
        ax[1].set_ylabel('tsne_qoi_2')
        ax[2].set_xlabel('tsne_1')
        ax[2].set_ylabel('tsne_2')

        fig.tight_layout()
        plt.show()


        pass
    def do_cluster_analysis(self, n_components):
        names_ = self.names
        data_ = self.data.df[names_]
        gmm = GaussianMixture(
                n_components=n_components,
                covariance_type='full',
                random_state=0).fit(data_)
        self.data.df['cluster_id'] = gmm.predict(data_)
        self.cluster_ids = list(set(self.data.df['cluster_id']))
        self.cluster_ids.sort()

        self.clusters = OrderedDict()
        for cluster_id in self.cluster_ids:
            self.clusters[cluster_id] = OrderedDict([
                ('cluster_id', cluster_id),
                ('N', self.data.df.loc[self.data.df['cluster_id'] == cluster_id].shape[0])
            ])
        
        n_clusters = len(self.cluster_ids)
        for i in range(n_clusters):
            self.clusters[i]['weight'] = gmm.weights_[i]
            self.clusters[i]['mean'] = gmm.means_[i,:]
            self.clusters[i]['covariance'] = gmm.covariances_[i,:]

        for i in range(n_clusters):
            self.clusters[i]['parameters'] = self._do_parameter_cluster_analysis(i)
            self.clusters[i]['qois'] = self._do_qoi_cluster_analysis(i)
    
    def _do_parameter_cluster_analysis(self, cluster_id):
        assert isinstance(cluster_id, int)
        assert cluster_id in self.cluster_ids

        data_ = self.data.df.loc[self.data.df['cluster_id'] == cluster_id]
        analysis_dict = OrderedDict()
        analysis_dict['mean'] = data_[self.configuration.parameter_names].mean()
        analysis_dict['std'] = data_[self.configuration.parameter_names].std()
        return analysis_dict

    def _do_qoi_cluster_analysis(self, cluster_id):
        assert isinstance(cluster_id, int)
        assert cluster_id in self.cluster_ids

        data_ = self.data.df.loc[self.data.df['cluster_id'] == cluster_id]
        analysis_dict = OrderedDict()
        analysis_dict['mean'] = data_[self.configuration.qoi_names].mean()
        analysis_dict['std'] = data_[self.configuration.qoi_names].std()
        return analysis_dict

    @staticmethod
    def plot_ellipse(position, covariance, ax=None, **kwargs):
        from matplotlib.matches import Ellipse

        if ax is None:
            fig, ax = plt.subplots(1,1)

        if covariance.shape == (2,2):
            U, s, Vt = np.linalg.svd(covariance)
            angle = np.degrees(np.arctan2(U[1,0],U[0,0]))
            width, height = 2 * np.sqrt(s)
        else:
            angle = 0
            width, height = 2 * np.sqrt(covariance)

        # draw ellipse
        for nsig in range(1,4):
            ax.add_patch(Ellipse(position,nsig*width,nsig*height,angle,**kwargs))

    @staticmethod
    def plot(gmm_obj, X, labels=None, ax=None, dpi=1200, filename=None,xlims=None,ylims=None):
        cluster_id = gmmi_obj.fit(X).predict(X)
        plt.close('all')

        if ax is None:
            fig, ax = plt. subplots(1,1)

        if isinstance(X, np.ndarray):
            x = X[:,0]
            y = X[:,1]

        elif isinstance(X, pd.DataFrame):
            x = X[X.columns[0]]
            y = X[X.columns[1]]

        ax.scatter(x,y,c=cluster_id,s=1,cmap='viridis',zorder=2,label=[k+1 for k in cluster_id])

        w_factor = 0.2 / gmm.weights_.max()
        for pos, covar, w in zip(gmm.means_,gmm.covariances_,gmm.weights_):
            plot_ellipse(pos,covar,alpha=w*w_factor,ax=ax)

        if labels is None:
            ax.set_xlabel(X.columns[0])
            ax.set_ylabel(X.columns[1])
        else:
            ax.set_xlabel(labels[0])
            ax.set_ylabel(labels[1])

        if xlims is not None:
            ax.set_xlim(xlims)
        if ylims is not None:
            ax.set_ylim(ylims)

        ax.legend()
        #ax.set(adjustable='box', aspect='equal')
        if filename is None:
            plt.show()
        else:
            fig.set_size_inches(5,5)
            fig.tight_layout()
            fig.savefig(filename,dpi=dpi)


