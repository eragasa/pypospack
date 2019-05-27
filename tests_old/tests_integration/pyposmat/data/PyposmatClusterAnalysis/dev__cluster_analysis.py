import os,time
from collections import OrderedDict

import numpy as np
import pandas as pd

from sklearn import preprocessing
from sklearn import manifold
from sklearn import neighbors
from sklearn import cluster
from sklearn import metrics

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import Imputer

from pypospack.exception import BadPreprocessorTypeException
from pypospack.exception import BadManifoldTypeException
from pypospack.exception import BadNearestNeighborTypeException
from pypospack.exception import BadClusterTypeException
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data.cluster_analysis import PyposmatClusterAnalysis
from pypospack.pyposmat.data.cluster_analysis import ecdf


import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

if __name__ == "__main__":

    pyposmat_src_dir = os.path.join(
        '..','..','..','..',
        'data_test','Ni__eam__born_exp_fs_00','data__Ni__eam__born_exp_fs_03'
    )
    pyposmat_configuration_fn = os.path.join(
        pyposmat_src_dir,'pyposmat.config.in'
    )
    pyposmat_data_fn = os.path.join(
        pyposmat_src_dir,'pyposmat.kde.5.out'
    )

    import numpy as np

    d = OrderedDict()
    d['configuration_fn'] = pyposmat_configuration_fn
    d['data_fn'] = pyposmat_data_fn
    d['include_parameters'] = True
    d['include_qois'] = True
    d['include_errors'] = False

    d['preprocessing'] = OrderedDict()
    d['preprocessing']['type'] = 'standard_scaler'
    d['preprocessing']['args'] = OrderedDict()
    d['preprocessing']['args']['copy'] = True
    d['preprocessing']['args']['with_mean'] = True
    d['preprocessing']['args']['with_std'] = True

    d['manifold'] = OrderedDict()
    d['manifold']['type'] = 'tsne'
    d['manifold']['args'] = OrderedDict()
    d['manifold']['args']['n_components'] = 2
    d['manifold']['args']['perplexity'] = 30
    d['manifold']['args']['early_exaggeration'] = 12
    d['manifold']['args']['learning_rate'] = 200
    d['manifold']['args']['n_iter'] = 5000
    d['manifold']['args']['n_iter_without_progress']=300,
    d['manifold']['args']['min_grad_norm'] =1e-7,
    #d['manifold']['args']['metric']='euclidean',
    d['manifold']['args']['init'] ='pca',
    d['manifold']['args']['verbose']=0,
    d['manifold']['args']['random_state']=None
    #method='barnes_hut'
    #angle=0.5

    d['neighbors'] = OrderedDict()
    d['neighbors']['type'] = 'ball_tree'
    d['neighbors']['kNN'] = 4
    d['neighbors']['args'] = OrderedDict()
    d['neighbors']['args']['leaf_size'] = 40
    d['neighbors']['args']['metric'] = 'minkowski'

    d['cluster'] = OrderedDict()
    d['cluster']['type'] = 'dbscan'
    d['cluster']['args'] = OrderedDict()
    d['cluster']['args']['eps'] = OrderedDict()
    d['cluster']['args']['eps']['NN'] = 3
    d['cluster']['args']['eps']['percentile'] = .99
    d['cluster']['args']['min_samples'] = 10
    d['cluster']['args']['metric'] = 'euclidean'
    d['cluster']['args']['metric_params'] = None
    d['cluster']['args']['algorithm'] = 'auto'
    d['cluster']['args']['leaf_size'] = 30
    d['cluster']['args']['p'] = None

    #o = PyposmatClusterAnalysis()
    o = PyposmatClusterAnalysis.init_from_ordered_dict(d)

    rd = OrderedDict()
    rd['include_parameters'] = d['include_parameters']
    rd['include_qois'] = d['include_qois']
    rd['include_errors'] = d['include_errors']
    rd['cluster'] = OrderedDict()
    rd['cluster']['names'] = o.cluster_names


    o.preprocess_data(d)
    rd['preprocessing'] = OrderedDict()
    rd['preprocessing']['type'] = d['preprocessing']['type']
    rd['preprocessing']['args'] = d['preprocessing']['args']
    rd['preprocessing']['names'] = o.scaled_names
    rd['preprocessing']['results'] = OrderedDict()
    #rd['preprocessing']['results']['mean'] = o.preprocessor.mean_
    #rd['preprocessing']['results']['sd'] = o.preprocessor.var_
    for i,v in enumerate(o.cluster_names):
        rd['preprocessing']['results'][v] = OrderedDict(
            [
                ('mean',o.preprocessor.mean_[i]),
                ('var',o.preprocessor.var_[i])
            ]
        )
    rd['preprocessing']['performance'] = OrderedDict()
    rd['preprocessing']['performance']['cpu_time'] = \
        o.cpu_time_preprocess

    if rd['preprocessing']['type'] is 'standard_scaler':
        for k,v in rd['preprocessing']['results'].items():
            print("{},{},{}".format(k,v['mean'],v['var']))

    o.calculate_manifold(d)

    rd['manifold'] = OrderedDict()
    rd['manifold']['type'] = d['manifold']['type']
    rd['manifold']['args'] = d['manifold']['args']
    rd['manifold']['results'] = OrderedDict()
    rd['manifold']['names'] = o.manifold_names
    rd['manifold']['performance'] = OrderedDict()
    rd['manifold']['performance']['cpu_time'] = \
        o.cpu_time_manifold


    o.calculate_kNN_analysis(d)
    o.calculate_clusters(d)

### PLOTTING RESULTS ----------------------------------------------------------

    fig1, (ax11,ax12) = plt.subplots(1,2)
    (nr,nc) = o.data.df.shape
    qe1,pe1 = ecdf(o.data.df['NN_1'])
    qe2,pe2 = ecdf(o.data.df['NN_2'])
    qe3,pe3 = ecdf(o.data.df['NN_3'])
    ax11.plot(qe1,pe1,lw=2,label='NN_1')
    ax11.plot(qe2,pe2,lw=2,label='NN_2')
    ax11.plot(qe3,pe3,lw=2,label='NN_3')
    ax11.set_xlabel('Quantile')
    ax11.set_ylabel('Cumulative probability')
    ax11.legend(fancybox=True, loc='right')
    ax12.hist(o.data.df['kNN_radius'])
    plt.show()

    fig2, (ax21,ax22,ax23) = plt.subplots(1,3)
    ax21.scatter(
        o.data.df[o.manifold_names[0]],
        o.data.df[o.manifold_names[1]],
        marker='.',
        s=1
    )

    cluster_colors = cm.spectral(
        o.data.df[o.data.df['cluster_id'] != -1]['cluster_id'].astype(float)\
            /o.n_clusters
    )

    ax22.scatter(
        o.data.df[o.data.df['cluster_id'] != -1][o.manifold_names[0]],
        o.data.df[o.data.df['cluster_id'] != -1][o.manifold_names[1]],
        marker='.',
        c=cluster_colors,
        s=1
    )
    ax23.scatter(
        o.data.df[o.data.df['cluster_id'] == -1][o.manifold_names[0]],
        o.data.df[o.data.df['cluster_id'] == -1][o.manifold_names[1]],
        marker='.',
        s=1
    )
    plt.show()

    from pandas.plotting import parallel_coordinates
    fig3, ax3 = plt.subplots(1,1)
    cluster_ids = set(o.data.df['cluster_id'])
    print(o.data.df.columns)
    for cluster_id in cluster_ids:
        parallel_coordinates(
            o.data.df[o.data.df['cluster_id'] == cluster_id][o.n_error_names],
            'Ni_fcc.E_coh.err.nerr'
        )
    plt.gca().legend_.remove()
    plt.show()
