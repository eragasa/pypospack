import os
from copy import deepcopy
from sklearn.mixture import GaussianMixture
import numpy as np
import pandas as pd

from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

class GmmAnalysis(object):

    def __init__(self,
                 configuration,
                 data,
                 names=None,
                 output_path='gmm_analysis',
                 max_components=20):
        self._initialize_configuration(configuration=configuration)
        self._initialize_data(data=data)
        self._initialize_names(names=names)
        self.names = deepcopy(names)
        self.output_path = output_path

        assert isinstance(max_components, int)
        self.max_components = max_components

        self.models = None
        self.aic_criteria = None
        self.bic_criteria = None
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
        else:
            raise TypeError

    def make_gmm_models(self, max_components=None):
        names_ = self.names
        data_ = self.data.df[names_]

        if max_components is not None:
            assert isinstance(max_components, int)
            self.max_components = max_components
        max_components_ = self.max_components

        self.models = {}
        n_components = [int(k) for k in np.arange(1, max_components_)]
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

        self.clusters = {}
        for cluster_id in self.cluster_ids:
            self.clusters[cluster_id] = {
                'cluster_id':cluster_id,
                'N':self.data.df.loc[self.data.df['cluster_id'] == cluster_id].shape[0]
            }

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
            print(pos)
            print(covar)
            print(w)
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


