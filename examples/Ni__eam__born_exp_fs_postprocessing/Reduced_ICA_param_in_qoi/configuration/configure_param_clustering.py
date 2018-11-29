from collections import OrderedDict
from pypospack.pyposmat.data.pipeline import PyposmatPipeline

pipeline_configuration = OrderedDict()

# define first segment (normalization)
pipeline_configuration[0] = OrderedDict() # int keys indicate step number
pipeline_configuration[0]['segment_type'] = 'preprocess'
pipeline_configuration[0]['function_calls'] = OrderedDict()
pipeline_configuration[0]['function_calls'][0]= OrderedDict()  # int keys allow multiple calls to same function
pipeline_configuration[0]['function_calls'][0]['function'] = 'normalize_standard_scaler'
pipeline_configuration[0]['function_calls'][0]['args'] = OrderedDict()
pipeline_configuration[0]['function_calls'][0]['args']['cols'] = ['param']
pipeline_configuration[0]['function_calls'][0]['args']['kwargs'] = OrderedDict()
pipeline_configuration[0]['function_calls'][0]['args']['kwargs']['standard_scaler'] = OrderedDict()
pipeline_configuration[0]['function_calls'][0]['args']['kwargs']['standard_scaler']['copy'] = False
pipeline_configuration[0]['function_calls'][0]['args']['kwargs']['standard_scaler']['with_mean'] = True
pipeline_configuration[0]['function_calls'][0]['args']['kwargs']['standard_scaler']['with_std'] = True

# define second segment (ICA transformation)
pipeline_configuration[1] = OrderedDict()
pipeline_configuration[1]['segment_type'] = 'pca'
pipeline_configuration[1]['function_calls'] = OrderedDict()
pipeline_configuration[1]['function_calls'][0]= OrderedDict()
pipeline_configuration[1]['function_calls'][0]['function'] = 'transform_ica'
pipeline_configuration[1]['function_calls'][0]['args'] = OrderedDict()
pipeline_configuration[1]['function_calls'][0]['args']['cols'] = ['n_param']
pipeline_configuration[1]['function_calls'][0]['args']['clusters'] = None
pipeline_configuration[1]['function_calls'][0]['args']['kwargs'] = OrderedDict()
pipeline_configuration[1]['function_calls'][0]['args']['kwargs']['ica'] = OrderedDict()
pipeline_configuration[1]['function_calls'][0]['args']['kwargs']['ica']['max_iter'] = 5000

# define third segment (clustering)
pipeline_configuration[2] = OrderedDict()
pipeline_configuration[2]['segment_type'] = 'cluster'
pipeline_configuration[2]['function_calls'] = OrderedDict()
pipeline_configuration[2]['function_calls'][0]= OrderedDict()
pipeline_configuration[2]['function_calls'][0]['function'] = 'cluster_kmeans'
pipeline_configuration[2]['function_calls'][0]['args'] = OrderedDict()
# throw out ica 0
pipeline_configuration[2]['function_calls'][0]['args']['abs_cols'] = ['ica_1', 'ica_2', 'ica_3', 'ica_4', 'ica_5', 'ica_6']
pipeline_configuration[2]['function_calls'][0]['args']['kwargs'] = OrderedDict()
pipeline_configuration[2]['function_calls'][0]['args']['kwargs']['kmeans'] = OrderedDict()
pipeline_configuration[2]['function_calls'][0]['args']['kwargs']['kmeans']['n_clusters'] = 10

# define fourth segment (plotting)
pipeline_configuration[3] = OrderedDict()
pipeline_configuration[3]['segment_type'] = 'plot'
pipeline_configuration[3]['function_calls'] = OrderedDict()
pipeline_configuration[3]['function_calls'][0]= OrderedDict()
pipeline_configuration[3]['function_calls'][0]['function'] = 'plot_by_cluster'
pipeline_configuration[3]['function_calls'][0]['args'] = OrderedDict()
pipeline_configuration[3]['function_calls'][0]['args']['x_axis'] = 'ica_1'
pipeline_configuration[3]['function_calls'][0]['args']['y_axis'] = 'ica_2'
pipeline_configuration[3]['function_calls'][0]['args']['filename'] = 'param_clusters_in_ica_space.png'


if __name__ == "__main__":
    pipeline = PyposmatPipeline()
    fn = __file__.replace('.py', '.in')
    pipeline.write_configuration(filename=fn,
                                 d=pipeline_configuration)
