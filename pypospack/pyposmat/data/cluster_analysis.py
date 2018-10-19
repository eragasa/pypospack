import time
from collections import OrderedDict

import numpy as np
import pandas as pd
import scipy as sp

from sklearn import preprocessing
from sklearn import manifold
from sklearn import neighbors
from sklearn import cluster
from sklearn import metrics

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import Imputer
from sklearn.decomposition import PCA

from pypospack.exceptions import BadPreprocessorTypeException
from pypospack.exceptions import BadManifoldTypeException
from pypospack.exceptions import BadNearestNeighborTypeException
from pypospack.exceptions import BadClusterTypeException
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.kde import Chiu1999_h


class PyposmatClusterSampler(object):

    def __init__(self):
        pass

class PyposmatPreprocessingPipeline(object):
    def __init__(self):
        pass

class PyposmatClusterAnalysis(object):
    def __init__(self, o_logger=None):
        if o_logger is None:
            raise NotImplementedError("backup logging not yet supported")
        else:
            self.log = o_logger
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
        if self.include_parameters is True:
            for v in self.parameter_names:
                scaled_names.append('{}.nparam'.format(v))
        if self.include_qois is True:
            for v in self.qoi_names:
                scaled_names.append('{}.nqoi'.format(v))
        if self.include_errors is True:
            for v in self.error_names:
                scaled_names.append('{}.nerr'.format(v))

        return scaled_names

    @property
    def n_parameter_names(self):
        n_parameter_names = [
            '{}.nparam'.format(v) for v in self.parameter_names
        ]
        return n_parameter_names

    @property
    def n_qoi_names(self):
        n_qoi_names = [
            '{}.nqoi'.format(v) for v in self.qoi_names
        ]
        return n_qoi_names

    @property
    def n_error_names(self):
        n_error_names = [
            # original error name is a.b.err, ignore .err
            '{}.nerr'.format(".".join(v.split('.')[:-1])) for v in self.error_names
        ]
        return n_error_names

    def init_from_json(json):
        pass

    def init_from_ordered_dict(d, o_logger=None):
        """
        constructor from a OrderedDict

        Arguments:
        ==========
        d(collections.OrderedDict):

        """

        o = PyposmatClusterAnalysis(o_logger=o_logger)

        o.configuration_fn = d['configuration_fn']
        try:
            o.read_configuration(filename=o.configuration_fn)
        except FileNotFoundError as e:
            print("ConfigurationError in PyposmatClusterAnalysis")
            for k,v in d.items():
                print("    {}={}".format(k,v))
            print("o.configuration_fn={}".format(
                o.configuration_fn)
            )
            raise e

        o.data_fn = d['data_fn']
        o.read_data(filename=o.data_fn)

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
        # quick patch to fix the duplicate axis issue
        if 'cluster_id' in list(o.data.df):
            o.data.df = o.data.df.drop(['cluster_id'], axis=1)
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
        self.nearest_neighbor_type = d['neighbors']['type']
        if self.nearest_neighbor_type == 'ball_tree':
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
                    qe[i],pe[i] = ecdf(self.data.df['NN_{}'.format(i)])
                kwargs['eps'] = np.interp(percentile,pe[kNN],qe[kNN])
            else:
                kwargs['eps'] = eps
        else:
            kwargs['eps'] = 0.5

        if 'min_samples' in d['cluster']['args']:
            min_samples = d['cluster']['args']['min_samples']

            self.data.df = pd.concat(
                [
                    self.data.df,
                    pd.DataFrame(
                        data = self.kNN_tree.query_radius(
                            self.data.df[self.manifold_names],
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
            self.data.df[self.manifold_names]
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
            self.df[
                self.parameter_names
                + self.error_names
                + self.qoi_names
            ]
        )

        self.data.df = pd.concat(
            [
                self.data.df,
                pd.DataFrame(
                    data=self.preprocessor.transform(
                        self.df[
                            self.parameter_names
                            + self.error_names
                            + self.qoi_names
                        ]
                    ),
                    columns=self.n_parameter_names\
                        + self.n_error_names\
                        + self.n_qoi_names
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
                self.data.df,
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

    def isValidPartition(self):
        # check the cholesky decomposition of each cluster
        # return True if no error
        # return False on linalg error
        for cluster_id in set(self.data.df['cluster_id']):
            cluster = self.data.df[self.data.df['cluster_id'] == cluster_id]
            cluster = cluster[self.parameter_names]
            try:
                _h = Chiu1999_h(cluster.values.T)
            except np.linalg.linalg.LinAlgError as e:
                self.log.write(e)
                self.log.write("cholesky failed on cluster {}".format(cluster_id))
                return False
            else:
                return True


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
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
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
    for cluster_id in cluster_ids:
        parallel_coordinates(
            o.data.df[o.data.df['cluster_id'] == cluster_id][o.n_error_names],
            'Ni_fcc.E_coh.nerr'
        )
    plt.gca().legend_.remove()
    plt.show()
