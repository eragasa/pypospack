import os
from collections import OrderedDict

import numpy as np
import pandas as pd
import pypospack.utils

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.manifold import TsneManifold

figure_path = os.path.join('gmm_analaysis_all',
                           'gmm_manifold.jpg')
figure_dpi = 1300
if __name__ == "__main__":
    # pypospack root directory
    pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()

    config_fn = os.path.join(
            pypospack_root_dir,
            'data','Si__sw__data','pareto_optimization_unconstrained',
            'pyposmat.config.in')
    data_fn = 'gmm_cluster.out'

    o_config = PyposmatConfigurationFile()
    o_config.read(filename=config_fn)
    
    o_data = PyposmatDataFile()
    o_data.read(filename=data_fn)

    for v in ['qois','free_parameters','all']:
        manifold = TsneManifold(pyposmat_configuration=config_fn,
                                pyposmat_data=data_fn)
        manifold.learn_manifold(v)

        if v == 'qois':
            o_data.df['TSNE_QOI_1'] = manifold.data.df['TSNE_1']
            o_data.df['TSNE_QOI_2'] = manifold.data.df['TSNE_2']
        elif v == 'free_parameters':
            o_data.df['TSNE_PARAM_1'] = manifold.data.df['TSNE_1']
            o_data.df['TSNE_PARAM_2'] = manifold.data.df['TSNE_2']
        elif v == 'all':
            o_data.df['TSNE_1'] = manifold.data.df['TSNE_1']
            o_data.df['TSNE_2'] = manifold.data.df['TSNE_2']
        else:
            raise ValueError()

    n_clusters = len(set(list(o_data.df['cluster_id'].values)))
    fig = plt.figure(figsize=(5,3))
    ax = []
    ax.append(fig.add_subplot(131))
    ax.append(fig.add_subplot(132))
    ax.append(fig.add_subplot(133))
    for i in range(len(ax)):
        ax[i].set_aspect('equal',adjustable='box')
        ax[i].xaxis.set_major_formatter(NullFormatter())
        ax[i].yaxis.set_major_formatter(NullFormatter())
    
    for i in range(n_clusters):
        ax[0].scatter(
                o_data.df['TSNE_QOI_1'].loc[o_data.df['cluster_id'] == i],
                o_data.df['TSNE_QOI_2'].loc[o_data.df['cluster_id'] == i],
                s=1)
        ax[1].scatter(
                o_data.df['TSNE_PARAM_1'].loc[o_data.df['cluster_id'] == i],
                o_data.df['TSNE_PARAM_2'].loc[o_data.df['cluster_id'] == i],
                s=1)
        ax[2].scatter(
                o_data.df['TSNE_1'].loc[o_data.df['cluster_id'] == i],
                o_data.df['TSNE_2'].loc[o_data.df['cluster_id'] == i],
                s=1)

    fig.tight_layout()
    fig.save_fig(figure_name,dpi=1200)
    plt.show()



