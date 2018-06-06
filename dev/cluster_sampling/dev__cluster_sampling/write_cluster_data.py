from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatDataFile
from mc_sampler_iterate_w_cluster import PyposmatClusterAnalysis

configuration_fn = 'data/pyposmat.config.in'
data_in_fn = 'data/pyposmat.kde.5.out'
data_out_fn = 'data/pyposmat.cluster.0.out'

d = OrderedDict()
d['configuration_fn'] = configuration_fn
d['data_fn'] = data_in_fn
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
d['manifold']['args']['n_iter_without_progress'] = 300,
d['manifold']['args']['min_grad_norm'] = 1e-7,
# d['manifold']['args']['metric']='euclidean',
d['manifold']['args']['init'] = 'pca',
d['manifold']['args']['verbose'] = 0,
d['manifold']['args']['random_state'] = None
# method='barnes_hut'
# angle=0.5

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


if __name__ == "__main__":
    obj_cluster_analysis = PyposmatClusterAnalysis.init_from_ordered_dict(d)
    print("type(obj_cluster_analysis)={}".format(str(type(obj_cluster_analysis))))
    obj_cluster_analysis.preprocess_data(d)
    
    obj_cluster_analysis.calculate_manifold(d)
    obj_cluster_analysis.calculate_kNN_analysis(d)
    obj_cluster_analysis.calculate_clusters(d)
    obj_cluster_analysis.data.write(data_out_fn)
