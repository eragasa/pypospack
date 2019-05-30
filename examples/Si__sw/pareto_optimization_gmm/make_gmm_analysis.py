"""
example on how to GMM

adapted from
Jake VanderPlass, Python Data Science Handbook
https://jakevdp.github.io/PythonDataScienceHandbook/05.12-gaussian-mixtures.html
"""
import os
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

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
    for i,(pos, covar, w) in enumerate(zip(gmm.means_,gmm.covariances_,gmm.weights_)):
        plot_ellipse(pos,covar,alpha=w*w_factor,ax=ax,fill=None)

    ax.set_xlabel(X.columns[0])
    ax.set_ylabel(X.columns[1])
    
    if filename is None:
        plt.show()
    else:
        fig.savefig(filename,dpi=dpi)

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
    gmm_output_dir = 'gmm_analysis'

    o_config = PyposmatConfigurationFile()
    o_config.read(filename=config_fn)

    o_data = PyposmatDataFile()
    o_data.read(filename=data_fn)
    o_data.create_normalized_errors(normalize_type='by_qoi_target',qoi_targets=o_config.qoi_targets)
    o_data.df['score'] = o_data.df[o_config.normalized_error_names].abs().sum(axis=1)

    name_1 = o_config.qoi_names[0]
    name_2 = o_config.qoi_names[1]
    data = o_data.df[[name_1,name_2]] 
    max_components = 21
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
    plot_fn='test_aic_bic.eps'
    plot_dpi=1200
    fig.savefig(plot_fn,dpi=plot_dpi)

    plt.close()

    filename = 'gmm_analysis.eps'
    plot_gmm(
            models[bic_n_components],
            data,
            filename=filename)

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

    parallel_plot_fn = 'gmm_parallel_plot.eps'
    o_plot.save_figure(filename=parallel_plot_fn)
    exit()
