""" wrapper for gaussian mixture models of sklearn"""

# author: Eugene J. Ragasa <eragasa@ufl.edu>
# License: MIT
import os

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from sklearn.mixture import GaussianMixture

import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

def plot_ellipse(position,covariance,ax=None,**kwargs):
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

def plot_gmm(gmm,X,label=True,ax=None,dpi=1200,filename=None):
    if ax is None:
        fig, ax = plt.subplots(1,1)

    cluster_id = gmm.fit(X).predict(X)

    if isinstance(X,np.ndarray):
        x = X[:,0]
        y = X[:,1]
    elif isinstance(X,pd.DataFrame):
        x = X[X.columns[0]]
        y = X[X.columns[1]]

    if label:
        ax.scatter(x,y,c=cluster_id,s=1,cmap='viridis',zorder=2)
    else:
        ax.scatter(x,y,s=40,zorder=2)

    w_factor = 0.2 / gmm.weights_.max()
    for pos, covar, w in zip(gmm.means_,gmm.covariances_,gmm.weights_):
        plot_ellipse(pos,covar,alpha=w*w_factor,ax=ax)

    ax.set_xlabel(X.columns[0])
    ax.set_ylabel(X.columns[1])

    if filename is None:
        plt.show()
    else:
        fig.savefig(filename,dpi=dpi)

class AbstractAnalysisClass():
    def __init__(self,
                 pyposmat_configuration,
                 pyposmat_data,
                 names):
        self.configuration = None
        self.data = None
        self.names = None
        self.initialize_configuration(configuration=pyposmat_configuration)
        self.initialize_data(data=pyposmat_data)
        self.initialize_names(names=names)
    
    def initialize_configuration(self,configuration):
        if isinstance(configuration,PyposmatConfigurationFile):
            self.configuration = configuration
        elif isinstance(configuration, str):
            self.configuration = PyposmatConfigurationFile()
            self.configuration.read(filename=configuration)
        else:
            msg = ("configuration must either be a path to a configuration "
                   "file or an instance of PyposmatConfigurationFile")
            raise TypeError(msg)

    def initialize_data(self, data):
        if isinstance(data, PyposmatDataFile):
            self.data = data
        elif isinstance(data, str):
            self.data = PyposmatDataFile()
            self.data.read(filename=data)
        else:
            msg = ("data must either be a path to a data file or an instance "
                    "of PyposmatDataFile")
            raise TypeError(msg)

    def initialize_names(self,names):
        if isinstance(names,list):
            self.names = names
        elif isinstance(names,str):
            if names == 'all':
                self.names = self.configuration.free_parameter_names\
                        + self.configuration.qoi_names
            elif names == 'free_parameters':
                self.names = self.configuration.free_parameter_names
            elif names == 'qois':
                self.names = self.configuration.qoi_names
            else:
                msg = ("valid names arguments are all,free_parameters,qois or "
                       "a list of names")
                raise ValueError(msg)
        else:
            msg = ("valid names arguments are all,free_parameters,qois or "
                   "a list of names")
            raise TypeError(msg)

class GmmInformationCriteriaAnalysis(AbstractAnalysisClass):

    def __init__(self,
                 pyposmat_configuration,
                 pyposmat_data,
                 names,
                 max_components=20):

        AbstractAnalysisClass.__init__(
                pyposmat_configuration=pyposmat_configuration,
                pyposmat_data=pyposmat-data,
                names=names)

        assert isinstance(max_components,int)
        assert max_components > 0
        self.max_components = max_components
        
        self.information_criteria = None

    def do_analysis(self):
        n_components = np.arange(1,max_components)
        self.models = [
                GaussianMixture(n_components=n,covariance_type='full',
                                random_state=0).fit(
                                    self.data.df[self.names]) \
                for n in n_components
                ]

        aic_criteria = [m.aic(self.data.df[self.names]) for m in self.models]
        bic_criteria = [m.bic(self.data.df[self.names]) for m in self.models]
        self.information_criteria = pd.DataFrame(data=OrderedDict([
                ('n_components',n_components),
                ('aic',aic_criteria)
                ('bic',bic_criteria)
                ]))

    def write(self,filename='aic_bic.csv'):
        with open(filename,'w') as f:
            f.write(self.information.criteria.to_csv)

    def plot(self,filename,dpi=1200):
        plt.close("all")

class GmmClusterAnalysis(AbstractAnalysisClass):

    def __init__(self,
                 pyposmat_configuration,
                 pyposmat_data,
                 names,
                 n_components):
        AbstractAnalysisClass.__init__(
                self,
                pyposmat_configuration=pyposmat_configuration,
                pyposmat_data=pyposmat_data,
                names=names)

        assert isinstance(n_components,int)
        assert n_components > 0
        self.n_components = n_components
        self.initialize_model(n_components=n_components)

    def initialize_model(self,n_components):
        self.model = GaussianMixture(n_components=n_components,
                                     covariance_type='full',
                                     random_state=0)
        self.model.fit(self.data.df[self.names])
        self.data.df['cluster_id'] = self.model.predict(self.data.df[self.names])
class GmmDataFile():
    pass

def gmm_aic_bic_analysis(
        config_fn,
        data_fn,
        names,
        output_directory='gmm_analysis',
        max_components=20):
    assert isinstance(config_fn,str)
    assert isinstance(data_fn,str)
    assert os.path.isfile(config_fn)
    assert os.path.isfile(data_fn)

    o_config = PyposmatConfigurationFile()
    o_config.read(filename=config_fn)

    o_data = PyposmatDataFile()
    o_data.read(filename=data_fn)
    o_data.create_normalized_errors(
            normalize_type='by_qoi_target',
            qoi_targets=o_config.qoi_targets)
    o_data.df['score'] = o_data.df[o_config.normalized_error_names].abs().sum(axis=1)

    data = o_data.df[names]

    n_components = np.arange(1,max_components)
    models = [GaussianMixture(n_components=n,covariance_type='full',random_state=0).fit(data) for n in n_components]

    # AIC analysis
    aic, aic_idx = min((val,idx) for (idx,val) in enumerate([m.aic(data) for m in models]))
    aic_n_components = n_components[aic_idx]
    aic_criteria = [m.aic(data) for m in models]
    # BIC analysis
    bic, bic_idx = min((val,idx) for (idx,val) in enumerate([m.bic(data) for m in models]))
    bic_n_components = n_components[bic_idx]
    bic_criteria = [m.bic(data) for m in models]

    #plot the criteria
    print('bic_n_components:{}'.format(bic_n_components))
    print('aic_n_components:{}'.format(aic_n_components))
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
    plot_fn=os.path.join(output_directory,'aic_bic_plot.eps')
    plot_gmm_aic_bic(
            filename=plot_fn,
            n_components=n_components,
            aic_criteria=aic_criteria,
            bic_criteria=bic_criteria,
            aic_n_components=aic_n_components,
            bic_n_components=bic_n_components)
    return aic_n_components, bic_n_components
