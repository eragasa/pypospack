"""
example on how to GMM

adapted from
Jake VanderPlass, Python Data Science Handbook
https://jakevdp.github.io/PythonDataScienceHandbook/05.12-gaussian-mixtures.html
"""
import os
import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
#import seaborn as sns; sns.set()
import numpy as np
import pandas as pd


from sklearn.mixture import GaussianMixture
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

def gmm_analysis(
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
    plot_fn=os.path.join(output_directory,'aic_bic_plot.eps')
    plot_gmm_aic_bic(
            filename=plot_fn,
            n_components=n_components,
            aic_criteria=aic_criteria,
            bic_criteria=bic_criteria,
            aic_n_components=aic_n_components,
            bic_n_components=bic_n_components)
    
def plot_gmm_aic_bic(
        filename,
        n_components,
        aic_criteria,
        bic_criteria,
        bic_n_components=None,
        aic_n_components=None,
        aic_color='blue',
        bic_color='red'):

    fig,ax = plt.subplots(1,1)
    
    plot_dpi=1200
    loc_legend = 'upper right'

    print('plot_fn:{}'.format(filename))
    #plot the AIC
    print('\tprinting the AIC')
    ax.plot(n_components,aic_criteria, label='AIC criteria',linestyle='-',color=aic_color)
    if aic_n_components is not None:
        ax.axvline(aic_n_components,label='AIC estimate',linestyle=':',color=aic_color)

    #plot the BIC
    print('\tprinting the BIC')
    ax.plot(n_components,bic_criteria, label='BIC criteria',linestyle='-',color=bic_color)
    if bic_n_components is not None:
        ax.axvline(bic_n_components,label='BIC estimate',linestyle=':',color=bic_color)

    x_lim_min = n_components.min()
    x_lim_max = n_components.max()
    ax.set_xlim(x_lim_min,x_lim_max)
    ax.legend(loc=loc_legend)
    ax.set_xlabel('number of components')
    ax.set_ylabel('criteria')

    plt.tight_layout()
    fig.savefig(filename,dpi=plot_dpi)
    plt.close('all')
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
    ref_config_fn = os.path.join(
            pypospack_root_dir,
            'data','Si__sw__data','reference_potentials',
            'pyposmat.config.in')
    ref_data_fn = os.path.join(
            pypospack_root_dir,
            'data','Si__sw__data','reference_potentials',
            'pyposmat.kde.1.out')

    o_config = PyposmatConfigurationFile()
    o_config.read(filename=config_fn)

    o_data = PyposmatDataFile()
    o_data.read(filename=data_fn)
    o_data.create_normalized_errors(
            normalize_type='by_qoi_target',
           qoi_targets=o_config.qoi_targets)
    print(o_config.normalized_error_names)
    print(o_data.df.columns)
    o_data.df['score'] = o_data.df[o_config.normalized_error_names].abs().sum(axis=1)

    name_1 = o_config.qoi_names[0]
    name_2 = o_config.qoi_names[1]
    data = o_data.df[[name_1,name_2]] 
    max_components = 21
    
    gmm_analysis(
        config_fn=config_fn,
        data_fn=data_fn,
        names=[name_1,name_2],
        output_directory='gmm_analysis',
        max_components=20)

    exit()
    n_components = np.arange(1,max_components)
    models = [GaussianMixture(n_components=n,covariance_type='full',random_state=0).fit(data) 
            for n in n_components]

    bic, idx = min((val, idx) for (idx, val) in enumerate([m.bic(data) for m in models]))
    bic_n_components = n_components[idx]

    aic, idx = min((val, idx) for (idx, val) in enumerate([m.aic(data) for m in models]))
    aic_n_components = n_components[idx]

    print("bic:{:.4f} @ {} components".format(bic,bic_n_components))
    print("aic:{:.4f} @ {} components".format(aic,aic_n_components))
    plt.close()
    fig, ax = plt.subplots(1,1)
    ax.plot(n_components, [m.bic(data) for m in models], label='BIC')
    ax.plot(n_components, [m.aic(data) for m in models], label='AIC')
    ax.legend(loc='best')
    ax.set_xlabel('n_components')
    fig.savefig(os.path.join('gmm_analysis','aic_bic_plt.eps'),dpi=1200)

    plt.close()

    filename = os.path.join('gmm_analysis','gmm_analysis.eps')
    plot_gmm(models[bic_n_components],data,filename=filename)

    from pypospack.pyposmat.visualization import PyposmatParallelCoordinatesPlot

    gmm = GaussianMixture(
            n_components=bic_n_components,
            covariance_type='full',
            random_state=0
            ).fit(data)
    o_data.df['cluster_id'] = gmm.predict(data)
    cluster_ids = set(o_data.df['cluster_id'])
    cluster_ids = list(cluster_ids)
    cluster_ids.sort()

    for cluster_id in cluster_ids:
        print("{},{}".format(
            cluster_id,
            o_data.df.loc[o_data.df['cluster_id'] == cluster_id].shape[0]))

    from matplotlib import cm
    n_clusters = len(cluster_ids)
    o_ref_data = PyposmatDataFile()
    o_ref_data.read(filename=ref_data_fn)
    n_reference_potentials =o_ref_data.df.shape[0]
    cm_subsection = np.linspace(0.0,1.0, n_clusters+n_reference_potentials)
    colors = [cm.jet(i) for i in cm_subsection]
    o_plot = PyposmatParallelCoordinatesPlot()
    for i, cluster_id in enumerate(cluster_ids):
        df = o_data.df.loc[o_data.df['cluster_id'] == cluster_id]
        o_plot.plot(
                config=o_config,
                data=o_data,
                label=cluster_id,
                nsmallest=20,
                linewidth=1,
                alpha=0.7,
                color=colors[i],
                cluster_id=cluster_id)

    # plot the reference potentials for comparison
    o_plot.plot_reference_potentials(
            config=ref_config_fn,
            data=ref_data_fn,
            linewidth=5)

    o_plot.save_figure(filename=os.path.join('gmm_analysis','gmm_parallel_plot.eps')) 
    exit()
