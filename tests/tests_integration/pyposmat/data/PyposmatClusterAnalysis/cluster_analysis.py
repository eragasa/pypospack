import time
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

class PyposmatClusterSampler(object):

    def __init__(self):
        pass

class PyposmatPreprocessingPipeline(object):
    def __init__(self):
        pass

class PyposmatClusterAnalysis(object):
    def __init__(self):
        self.configuration_fn = None
        self.data_fn = None

        self.data = None
        self.configuration = None

        self.parameter_names = None
        self.qoi_names = None
        self.error_names = None

        self.include_parameters = True
        self.include_qois = False
        self.include_errors = False

        self.manifold_names = None

    @property
    def cluster_names(self):
        # determine which variables are included in the analysis
        cluster_names = []
        if self.include_parameters:
            cluster_names += self.parameter_names
        if self.include_qois:
            cluster_names += self.qoi_names
        if self.include_errors:
            cluster_names += self.error_names

        return cluster_names

    @property
    def scaled_names(self):

        scaled_names = []
        if o.include_parameters is True:
            for v in o.parameter_names:
                scaled_names.append('{}.nparam'.format(v))
        if o.include_qois is True:
            for v in o.qoi_names:
                scaled_names.append('{}.nqoi'.format(v))
        if o.include_errors is True:
            for v in o.error_names:
                scaled_names.append('{}.nerr'.format(v))

        return scaled_names

    def init_from_json(json):
        pass

    def init_from_ordered_dict(d):
        """
        constructor from a OrderedDict

        Arguments:
        ==========
        d(collections.OrderedDict):

        """

        o = PyposmatClusterAnalysis()

        o.configuration_fn = d['configuration_fn']
        o.read_configuration(filename=o.configuration_fn)

        o.data_fn = d['data_fn']
        o.read_data(filename=o.data_fn)

        #
        if 'include_parameters' in d:
            o.include_parameters = d['include_parameters']
        else:
            d['include_parameters'] = o.include_parameters
        if 'include_qois' in d:
            o.include_qois = d['include_qois']
        else:
            d['include_qois'] = o.include_qois
        if 'include_errors' in d:
            o.include_errors = d['include_errors']
        else:
            d['include_errors'] = o.include_errors

        if 'normalizer' in d:
            o.normalize_data(**d['normalizer'])

        if 'normalizer_type' in d:
            o.normalizer_type = d['normalizer_type']
            o.normalize_data(o.normalizer_type)

        return o

    def to_dict(self):
        d = OrderedDict()
        d['configuration_fn'] = o.configuration_fn
        d['data_fn'] = o.data_fn
        d['rescaler_type'] = o.rescaler_type
        d['analysis_on']['include_parameters'] = o.include_parameters
        d['analysis_on']['include_errors'] = o.include_errors

    def to_json(self):
        pass

    def to_ordered_dict(self):
        d = OrderedDict()
        d['pyposmat_config_fn'] = self.configuration_fn
        d['pyposmat_data_fn'] =self.data_fn

    def read_configuration(self,filename):
        """
        read in pyposmat configuration file

        Argum ents:
        ==========
        filename(str):
        """
        self.configuration_fn = filename
        self.configuration = PyposmatConfigurationFile()
        self.configuration.read(filename=filename)

    def read_data(self,filename):
        """
        read in pyposmat data filename

        Arguments:
        ==========
        filename(str):
        """

        self.data_fn = filename
        self.data = PyposmatDataFile()
        self.data.read(filename)

        self.parameter_names = self.data.parameter_names
        self.qoi_names = self.data.qoi_names
        self.error_names = self.data.error_names
        self.df = self.data.df

    def preprocess_data(self,d):
        self.preprocessor_type = d['preprocessing']['type']
        self.preprocessor = None

        time_start = time.time()
        if self.preprocessor_type == 'standard_scaler':
            self._preprocess_data__standard_scaler(d)
        else:
            raise BadPreprocessorTypeException()
        time_end = time.time()
        self.cpu_time_preprocess = time_end-time_start

    def calculate_manifold(self,d):
        time_start = time.time()

        self.manifold = None
        self.manifold_type = d['manifold']['type']
        if self.manifold_type == 'tsne':
            self._manifold_tsne(d)
        else:
            raise BadManifoldTypeException()

        time_stop = time.time()
        self.cpu_time_manifold = time_stop - time_start

    def calculate_kNN_analysis(self,d):
        start_time = time.time()
        self.kNN_tree = None
        self.kNN_names = None
        self.kNN = d['neighbors']['kNN']
        o.nearest_neighbor_type = d['neighbors']['type']
        if o.nearest_neighbor_type == 'ball_tree':
            self._kNN__ball_tree(d)
        else:
            raise BadNearestNeighborTypeException()
        stop_time = time.time()
        self.cpu_time_kNN = stop_time - start_time


    def calculate_clusters(self,d):

        time_start = time.time()
        self.clusterer = None
        self.cluster_type = d['cluster']['type']

        if self.cluster_type == 'dbscan':
            self._cluster_dbscan(d)
        else:
            raise BadClusterTypeException()
        time_stop = time.time()
        self.cpu_time_cluster = time_stop - time_start


    def _cluster_dbscan(self,d):
        kwargs = OrderedDict()
        if 'eps' in d['cluster']['args']:
            eps = d['cluster']['args']['eps']
            if isinstance(eps,dict):
                kNN = eps['NN']
                percentile = eps['percentile']

                qe = OrderedDict()
                pe = OrderedDict()
                for i in range(self.kNN):
                    qe[i],pe[i] = ecdf(o.data.df['NN_{}'.format(i)])
                kwargs['eps'] = np.interp(percentile,pe[kNN],qe[kNN])
            else:
                kwargs['eps'] = eps
        else:
            kwargs['eps'] = 0.5

        if 'min_samples' in d['cluster']['args']:
            min_samples = d['cluster']['args']['min_samples']

            o.data.df = pd.concat(
                [
                    o.data.df,
                    pd.DataFrame(
                        data = self.kNN_tree.query_radius(
                            o.data.df[self.manifold_names],
                            r = kwargs['eps'],
                            count_only = True
                        ),
                        columns = ['kNN_radius']
                    )
                ],
                axis = 1
            )
            kwargs['min_samples'] = min_samples

        self.clusterer = cluster.DBSCAN(**kwargs)
        data = self.clusterer.fit(
            o.data.df[o.manifold_names]
        )

        cluster_labels = data.labels_

        self.n_clusters = len(set(cluster_labels))
        self.n_clusters = self.n_clusters - (1 if -1 in cluster_labels else 0)

        self.data.df = pd.concat(
            [
                self.data.df,
                pd.DataFrame(
                    data=list(cluster_labels),
                    columns=['cluster_id']
                )
            ],
            axis = 1
        )

    def _preprocess_data__standard_scaler(self,d):
        # default keyword arguments
        kwargs = OrderedDict()
        kwargs['copy'] = True
        kwargs['with_mean'] = True
        kwargs['with_std'] = True

        # change default arguments
        for k,v in d['preprocessing']['args'].items():
            kwargs[k] = v

        # add the default arguments to the
        for k,v in kwargs.items():
            d['preprocessing']['args'][k] = v

        #pipeline_preprocessor = preprocessing.StandardScaler()
        self.preprocessor = preprocessing.StandardScaler(**kwargs)
        self.preprocessor.fit(
            self.df[self.cluster_names]
        )

        self.data.df = pd.concat(
            [
                self.data.df,
                pd.DataFrame(
                    data=self.preprocessor.transform(
                        self.df[self.cluster_names]
                    ),
                    columns=self.scaled_names
                )
            ],
            axis = 1
        )

    def _manifold_tsne(self,d):
        # default keyword arguments
        kwargs = OrderedDict()
        kwargs['n_components'] = 2
        kwargs['perplexity'] = 30
        kwargs['early_exaggeration'] = 12
        kwargs['learning_rate'] = 200
        kwargs['n_iter'] = 1000
        kwargs['n_iter_without_progress']=300,
        kwargs['min_grad_norm'] =1e-7,
        #kwargs['metric']='euclidean',
        kwargs['init'] ='pca',
        kwargs['verbose']=0,
        kwargs['random_state']=None

        # change default arguments
        for k,v in d['manifold']['args'].items():
            kwargs[k] = v

        # add the default arguments to the
        for k,v in kwargs.items():
            d['manifold']['args'][k] = v

        self.manifold_names = []
        for i in range(kwargs['n_components']):
            self.manifold_names.append('tsne_{}'.format(i))

        #TODO: need to fix
        self.manifold =  manifold.TSNE()
        #pipeline_manifold =  manifold.TSNE(**kwargs)

        self.data.df = pd.concat(
            [
                o.data.df,
                pd.DataFrame(
                    data = self.manifold.fit_transform(
                        self.data.df[self.scaled_names]
                    ),
                    columns = self.manifold_names
                )
            ],
            axis = 1
        )

    def _kNN__ball_tree(self,d):
        kwargs = OrderedDict()
        kwargs['leaf_size'] = 40
        kwargs['metric'] = 'minkowski'

        for k,v in d['neighbors']['args'].items():
            kwargs[k] = v

        for k,v in kwargs.items():
            d['neighbors']['args'][k] = v

        self.kNN_tree = neighbors.BallTree(
            self.data.df[self.manifold_names],
            **kwargs
        )

        self.kNN_names = ['NN_{}'.format(i) for i in range(self.kNN)]

        self.data.df = pd.concat(
            [
                self.data.df,
                pd.DataFrame(
                    data = self.kNN_tree.query(
                        X = self.data.df[self.manifold_names],
                        k = self.kNN,
                        return_distance = True
                    )[0],
                    columns = self.kNN_names
                )
            ],
            axis = 1
        )

def ecdf(sample):
    # convert sample to a numpy array, if it isn't already
    sample = np.atleast_1d(sample)

    # find the unique values and their corresponding counts
    quantiles, counts = np.unique(sample, return_counts=True)

    # take the cumulative sum of the counts and divide by the sample size to
    # get the cumulative probabilities between 0 and 1
    cumprob = np.cumsum(counts).astype(np.double) / sample.size

    return quantiles, cumprob

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

class PyposmatSilhouettePlot(object):
    """ Plot Silhouette Plots

    References:
    ===========
    http://scikit-learn.org/stable/auto_examples/cluster/plot_kmeans_silhouette_analysis.html#sphx-glr-auto-examples-cluster-plot-kmeans-silhouette-analysis-py
    """
    def __init__(self):
       self.fig = None
       self.silhouette_plot = None
       self.cluster_plot = None
       self.fig, (self.silhouette_plot,self.cluster_plot) = plt.subplots(1,2)

       # 1st subplot is silhouette plot_kmeans_silhouette_analysis

if __name__ == "__main__":
    import os

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
    d['include_qois'] = False
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


    #ax12.plot(o.data.df['NN_1'].sort_values(),range(nr),lw=2,label='NN_1')
    #ax12.plot(o.data.df['NN_2'].sort_values(),range(nr),lw=2,label='NN_2')
    #ax12.plot(o.data.df['NN_3'].sort_values(),range(nr),lw=2,label='NN_3')
    #ax12.set_xlabel('Quantile')
    #ax12.set_ylabel('Cumulative probability')
    #ax12.legend(fancybox=True, loc='right')


    print('scaling data')
    time_start = time.time()
    time_start_preprocessing = time.time()
    preprocessor_type = d['preprocessing']['type']
    pipeline_preprocessor = None

    if preprocessor_type == 'standard_scaler':

        # default keyword arguments
        kwargs = OrderedDict()
        kwargs['copy'] = True
        kwargs['with_mean'] = True
        kwargs['with_std'] = True

        # change default arguments
        for k,v in d['preprocessing']['args'].items():
            kwargs[k] = v

        # add the default arguments to the
        for k,v in kwargs.items():
            d['preprocessing']['args'][k] = v

        #pipeline_preprocessor = preprocessing.StandardScaler()
        pipeline_preprocessor = preprocessing.StandardScaler(**kwargs)
        pipeline_preprocessor.fit(o.df[o.cluster_names])

        o.data.df = pd.concat(
            [
                o.data.df,
                pd.DataFrame(
                    data=pipeline_preprocessor.transform(o.df[o.cluster_names]),
                    columns=o.scaled_names
                )
            ],
            axis = 1
        )
        print(o.data.df.shape)

        time_end_preprocessing = time.time()
    else:
        raise BadPreprocessorTypeException()


    # add manifold learner to pipeline
    print('learning_manifold')
    manifold_time_start = time.time()
    pipeline_manifold = None
    data_manifold = None
    manifold_names = None
    o.manifold_type = d['manifold']['type']
    if o.manifold_type == 'tsne':

        # default keyword arguments
        kwargs = OrderedDict()
        kwargs['n_components'] = 2
        kwargs['perplexity'] = 30
        kwargs['early_exaggeration'] = 12
        kwargs['learning_rate'] = 200
        kwargs['n_iter'] = 1000
        kwargs['n_iter_without_progress']=300,
        kwargs['min_grad_norm'] =1e-7,
        #kwargs['metric']='euclidean',
        kwargs['init'] ='pca',
        kwargs['verbose']=0,
        kwargs['random_state']=None

        # change default arguments
        for k,v in d['manifold']['args'].items():
            kwargs[k] = v

        # add the default arguments to the
        for k,v in kwargs.items():
            d['manifold']['args'][k] = v

        manifold_names = []
        for i in range(kwargs['n_components']):
            manifold_names.append('tsne_{}'.format(i))

        #TODO: need to fix
        pipeline_manifold =  manifold.TSNE()
        #pipeline_manifold =  manifold.TSNE(**kwargs)

        o.data.df = pd.concat(
            [
                o.data.df,
                pd.DataFrame(
                    data = pipeline_manifold.fit_transform(
                        o.data.df[o.scaled_names]
                    ),
                    columns = manifold_names
                )
            ],
            axis = 1
        )
        manifold_time_stop = time.time()

    else:
        raise BadManifoldTypeException()

    print(o.data.df.shape)
    for n in manifold_names:
        assert n in list(o.data.df.columns.values)
    rd['manifold'] = OrderedDict()
    rd['manifold']['type'] = d['manifold']['type']
    rd['manifold']['args'] = d['manifold']['args']
    rd['manifold']['results'] = OrderedDict()
    rd['manifold']['axes'] = list(manifold_names)
    rd['manifold']['performance'] = OrderedDict()
    rd['manifold']['performance']['cpu_time'] = \
        manifold_time_stop - manifold_time_start

    print('nearest neighbors analysis')
    kNN_start_time = time.time()
    kNN_tree = None
    kNN_names = None
    kNN = 4
    o.nearest_neighbor_type = d['neighbors']['type']
    if o.nearest_neighbor_type == 'ball_tree':
        kwargs = OrderedDict()
        kwargs['leaf_size'] = 40
        kwargs['metric'] = 'minkowski'

        for k,v in d['neighbors']['args'].items():
            kwargs[k] = v

        for k,v in kwargs.items():
            d['neighbors']['args'][k] = v

        kNN_tree = neighbors.BallTree(
            o.data.df[manifold_names],
            **kwargs
        )

        kNN_names = ['NN_{}'.format(i) for i in range(kNN)]

        o.data.df = pd.concat(
            [
                o.data.df,
                pd.DataFrame(
                    data = kNN_tree.query(
                        X = o.data.df[manifold_names],
                        k = kNN,
                        return_distance = True
                    )[0],
                    columns = kNN_names
                )
            ],
            axis = 1
        )

    else:
        raise BadNearestNeighborTypeException()

    print(o.data.df[kNN_names])
    print(o.data.df.shape)
    import matplotlib.pyplot as plt

    # Here we need to select an algorithm to determine the EPS based on the
    # distribution of the nearest neighbors.  We calculate the empirical
    # cumulative distribution function.

    # For now, we do visual inspection of the empirical CDF
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

    ax12.plot(o.data.df['NN_1'].sort_values(),range(nr),lw=2,label='NN_1')
    ax12.plot(o.data.df['NN_2'].sort_values(),range(nr),lw=2,label='NN_2')
    ax12.plot(o.data.df['NN_3'].sort_values(),range(nr),lw=2,label='NN_3')
    ax12.set_xlabel('Quantile')
    ax12.set_ylabel('Cumulative probability')
    ax12.legend(fancybox=True, loc='right')
    plt.show()

    d['cluster']['args']['eps'] = np.interp(.95,pe3,qe3)
    print('cluster_eps={}'.format(d['cluster']['args']['eps']))
    stop_time_NN = time.time()

     # Here we estimate the number minimum point cluster leaf_size
    o.data.df = pd.concat(
        [
            o.data.df,
            pd.DataFrame(
                data = kNN_tree.query_radius(
                    o.data.df[manifold_names],
                    r = d['cluster']['args']['eps'],
                    count_only = True
                ),
                columns = ['kNN_radius']
            )
        ],
        axis = 1
    )

    fig2, ax2 = plt.subplots(1,1)
    ax2.hist(o.data.df['kNN_radius'])
    plt.show()

    # add clustering
    print('learning clusters')
    time_start = time.time()
    pipeline_cluster = None
    data_cluster = None

    o.cluster_type = d['cluster']['type']
    if o.cluster_type == 'dbscan':
        kwargs = d['cluster']['args']
        pipeline_cluster = cluster.DBSCAN(**kwargs)
        data_cluster = pipeline_cluster.fit(o.data.df[manifold_names])
        print(data_cluster)
    else:
        raise BadClusterTypeException()
    time_stop = time.time()
    if 'results' not in d['cluster']:
        d['cluster']['results'] = OrderedDict()
    #core_samples_mask = np.zeros_like(data_cluster.labels_, dtype=bool)
    #core_samples_mask[data_cluster.core_samples_indices_] = True
    cluster_labels = data_cluster.labels_

    n_clusters = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
    print("n_clusters:{}".format(n_clusters))


    silhouette_avg = None
    silhouette_scores = None
    if n_clusters > 0:
        silhouette_avg = metrics.silhouette_score(
            o.data.df[manifold_names],
            cluster_labels
        )
        silhouette_scores = metrics.silhouette_samples(
            o.data.df[manifold_names],
            cluster_labels
        )
        print("silhouette_score_avg:{}".format(silhouette_scores))

        for i in range(n_clusters):
            print(
                "cluster_id={}".format(i)
            )
            print(
                "\tsize={}".format(
                    silhouette_scores[cluster_labels == i].size
                )
            )
            print(
                "\tsilhouette_score={}".format(
                    silhouette_scores[cluster_labels == i].max()
                )
            )

    import matplotlib.pyplot as plt
    cluster_colors = cm.spectral(cluster_labels.astype(float)/n_clusters)
    plt.scatter(
        o.data.df[manifold_names[0]].iloc[cluster_labels == -1],
        o.data.df[manifold_names[1]].iloc[cluster_labels == -1],
        marker='.',
        c=cluster_colors,
        s=1
    )
    plt.show()
        # cluster data
