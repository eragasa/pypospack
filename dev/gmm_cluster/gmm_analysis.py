import os, shutil
from copy import deepcopy
from collections import OrderedDict

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
from sklearn.mixture import GaussianMixture

from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

class BaseAnalysis(object):
    def __init__(self,
                 configuration,
                 data,
                 output_path=None):
        self.configuration = None
        self.data = None
        self.output_path = None

        self._initialize_configuration(configuration=configuration)
        self._initialize_data(data=data)
        self._initialize_output_path(path=output_path)

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

    def _initialize_output_path(self, path):
        if path is None:
            self.output_path = None
        elif isinstance(path, str):
            if os.path.isdir(path):
                shutil.rmtree(path)
            os.mkdir(path) 
            self.output_path = path
        else:
            raise TypeError

class GmmAnalysis(BaseAnalysis):

    def __init__(self,
                 configuration,
                 data,
                 names=None,
                 output_path='gmm_analysis',
                 max_components=20):
        BaseAnalysis.__init__(self,
                              configuration=configuration,
                              data=data,
                              output_path=output_path)
        
        self._initialize_names(names=names)

        assert isinstance(max_components, int)
        self.max_components = max_components

        self.models = None
        self.aic_criteria = None
        self.bic_criteria = None
        self.cluster_ids = None

        self.ic_plot_path = None

        self.manifold_type = None

    def _initialize_names(self, names):
        if isinstance(names, list):
            self.names = list(names)
        elif isinstance(names, str):
            if names == 'qois':
                self.names = self.configuration.qoi_names
            elif names == 'parameters':
                self.names = self.configuration.parameter_names
            elif names == 'all':
                parameter_names_ = self.configuration.parameter_names
                qoi_names_ = self.configuration.qoi_names
                self.names = parameter_names_ + qoi_names_
            else:
                m = 'valid values for names is qoi, parameters, or all'
                raise ValueError(m)
        else:
            raise TypeError

    def do_ic_analysis(self, max_components=None):
        if max_components is not None:
            self.max_components = max_components

        self.make_gmm_models()
        self.do_aic_analysis()
        self.do_bic_analysis()

        self.plot_ic_path = os.path.join(self.output_path, 'ic_plot.png')
        self.plot_ic(path=self.plot_ic_path)

    def plot_ic(self, path=None):
        x = [k+1 for k in range(self.max_components)]
        x_min = 1
        x_max = self.max_components
        aic = [v['aic'] for k, v in self.models.items()]
        bic = [v['bic'] for k, v in self.models.items()]
        
        plt.close()
        fig, ax = plt.subplots(1, 1)
        ax.plot(x, aic, label='aic')
        ax.plot(x, bic, label='bic')
        #ax.hline(aic_n_components)
        #ax.hline(bic_n_components)
        ax.legend()
        ax.set_xlim([x_min, x_max])
        if path is None:
            plt.show()
        else:
            plt.savefig(path, bbox_inches='tight')


    def make_gmm_models(self, max_components=None):
        names_ = self.names
        data_ = self.data.df[names_]

        if max_components is not None:
            assert isinstance(max_components, int)
            self.max_components = max_components
        max_components_ = self.max_components

        self.models = {}
        n_components = [k+1 for k in range(max_components_)]
        for n in n_components:
            self.models[n] = {}
            self.models[n]['obj'] = GaussianMixture(n_components=n,
                                                    covariance_type='full',
                                                    random_state=0).fit(data_)

    def do_aic_analysis(self):
        # AIC analysis
        models_ = self.models
        names_ = self.names
        data_ = self.data.df[names_]
       
        for k in self.models:
            aic = self.models[k]['obj'].aic(data_)
            self.models[k]['aic'] = aic
        aic, aic_n_components = min([(v['aic'], k) for k, v in self.models.items()])
        
        self.aic_criteria = {
            'min_components':int(aic_n_components),
            'min_value':float(aic)
        }

    def do_bic_analysis(self):
        models_ = self.models
        names_ = self.names
        data_ = self.data.df[names_]
        
        # BIC analysis
        for k in self.models:
            bic = self.models[k]['obj'].bic(data_)
            self.models[k]['bic'] = bic
        bic, bic_n_components = min([(v['bic'], k) for k, v in self.models.items()])
        
        self.bic_criteria = {
            'min_components':int(bic_n_components),
            'min_value':float(bic)
        }

    def make_manifold(self, 
                      manifold_type='tsne', 
                      n_components=2, 
                      names='all'):

        if isinstance(names, str):
            if names == 'parameters':
                names_ = self.configuration.parameter_names
            elif names == 'qois':
                names_ = self.configuration.qoi_names
            elif names == 'all':
                parameter_names_ = self.configuration.parameter_names
                qoi_names_ = self.configuration.qoi_names
                names_ = parameter_names + qoi_names_
            else:
                raise ValueError
        else:
            raise TypeError

        data_ = self.data.df[names_]

        if manifold_type == 'tsne':
            self.manifold_type = manifold_type
            manifold = TSNE(n_components=n_components)
        else:
            raise ValueError
        manifold_components = manifold.fit_transform(data_)

        if names == 'parameters':
            for i in range(n_components):
                column_name = '{}_param_{}'.format(manifold_type, i+1)
                self.data.df[column_name] = manifold_components[:,i]
        elif names == 'qois':
            for i in range(n_components):
                column_name = '{}_param_{}'.format(manifold_type, i+1)
                self.data.df[column_name] = manifold_components[:,i]
        elif names == 'all':
            for i in range(n_components):
                column_name = '{}_{}'.format(manifold_type, i+1)
                self.data.df[column_name] = manifold_components[:,i]
        else:
            raise ValueError

    
    def plot_gmm_analysis(self, n_components):

        fig, ax = plt.subplots(1,3)

        print('learning parameter manifold...')
        self.make_manifold(names='parameters')
        print('learning qoi manifold...')
        self.make_manifold(names='qois')
        print('learning total manifold...')
        self.make_manifold(names='all')

        self.do_cluster_analysis(n_components=n_components)
        for i in self.cluster_ids:
            col_name_1 = '{}_param_{}'.format(self.manifold_type, 1)
            col_name_2 = '{}_param_{}'.format(self.manifold_type, 2)
            ax[0].scatter(
                self.data.df[col_name_1].loc[self.data.df['cluster_id']==i],
                self.data.df[col_name_2].loc[self.data.df['cluster_id']==i],
                s=1.)

            col_name_1 = '{}_qoi_{}'.format(self.manifold_type, 1)
            col_name_2 = '{}_qoi_{}'.format(self.manifold_type, 2)
            ax[1].scatter(
                self.data.df[col_name_1].loc[self.data.df['cluster_id']==i],
                self.data.df[col_name_2].loc[self.data.df['cluster_id']==i],
                s=1.)
            
            col_name_1 = '{}_{}'.format(self.manifold_type, 1)
            col_name_2 = '{}_{}'.format(self.manifold_type, 2)
            ax[2].scatter(
                self.data.df[col_name_1].loc[self.data.df['cluster_id']==i],
                self.data.df[col_name_2].loc[self.data.df['cluster_id']==i],
                s=1.)

        xlabel_str = '{}_param_{}'.format(self.manifold_type, 1)
        ylabel_str = '{}_param_{}'.format(self.manifold_type, 2)
        ax[0].set_xlabel(xlabel_str)
        ax[0].set_ylabel(ylabel_str)

        xlabel_str = '{}_qoi_{}'.format(self.manifold_type, 1)
        ylabel_str = '{}_qoi_{}'.format(self.manifold_type, 2)
        ax[1].set_xlabel(xlabel_str)
        ax[1].set_ylabel(ylabel_str)

        xlabel_str = '{}_{}'.format(self.manifold_type, 1)
        ylabel_str = '{}_{}'.format(self.manifold_type, 2)
        ax[2].set_xlabel(xlabel_str)
        ax[2].set_ylabel(ylabel_str)

        fig.tight_layout()
        plt.show()

    def do_cluster_analysis(self, n_components):
        names_ = self.names
        data_ = self.data.df[names_]
        self.gmm = GaussianMixture(
                n_components=n_components,
                covariance_type='full',
                random_state=0).fit(data_)

        self.data.df['cluster_id'] = self.gmm.predict(data_)
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
            self.clusters[i]['weight'] = self.gmm.weights_[i]
            self.clusters[i]['mean'] = self.gmm.means_[i,:]
            self.clusters[i]['covariance'] = self.gmm.covariances_[i,:]

        for i in range(n_clusters):
            self.clusters[i]['parameters'] = self._do_parameter_cluster_analysis(i)
            self.clusters[i]['qois'] = self._do_qoi_cluster_analysis(i)

    def table_cluster_qois(self):
        qoi_names_ = self.configuration.qoi_names

        header_row = ['cluster_id', ''] + [p for p in qoi_names_]
        
        mu_rows = []
        std_rows = []
        for k, gmm in self.clusters.items():
            mu_rows.append(
                    [k, 'mu'] + [gmm['qois']['mean'][p] for p in qoi_names_]
            )
            std_rows.append(
                    [k, 'sigma'] + [gmm['qois']['std'][p] for p in qoi_names_]
            )

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

        n_parameters = len(self.configuration.parameter_names)
        n_qois = len(self.configuration.qoi_names)
        n_ttl = n_parameters + n_qois
        qoi_map = range(n_parameters, n_ttl)
        data_ = self.data.df.loc[self.data.df['cluster_id'] == cluster_id]
        analysis_dict = OrderedDict()
        analysis_dict['mean'] = data_[self.configuration.qoi_names].mean()
        analysis_dict['std'] = data_[self.configuration.qoi_names].std()
        analysis_dict['mean_gmm'] = self.gmm.means_[cluster_id, qoi_map]
        analysis_dict['std_gmm'] = np.sqrt(np.array(
                [self.gmm.covariances_[cluster_id,k,k] for k in qoi_map]
        ))

        from numpy.linalg import  eig
        cov_matrix = self.gmm.covariances_[cluster_id]
        qoi_cov_matrix = cov_matrix[np.ix_(qoi_map, qoi_map)]
        eigval, eigvec = eig(qoi_cov_matrix)
        print(eigval, eigvec)
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


