from collections import OrderedDict
from pypospack.pyposmat.data.pipeline import PyposmatPipeline

pipeline_configuration = OrderedDict()


# define first segment (plotting)
pipeline_configuration[0] = OrderedDict()
pipeline_configuration[0]['segment_type'] = 'plot'
pipeline_configuration[0]['function_calls'] = OrderedDict()
pipeline_configuration[0]['function_calls'][0]= OrderedDict()
pipeline_configuration[0]['function_calls'][0]['function'] = 'plot_by_cluster'
pipeline_configuration[0]['function_calls'][0]['args'] = OrderedDict()
pipeline_configuration[0]['function_calls'][0]['args']['x_axis'] = 'pca_0'
pipeline_configuration[0]['function_calls'][0]['args']['y_axis'] = 'pca_1'
pipeline_configuration[0]['function_calls'][0]['args']['filename'] = 'err_clusters_in_param_pca_space.png'

if __name__ == "__main__":
    pipeline = PyposmatPipeline()
    fn = __file__.replace('.py', '.in')
    pipeline.write_configuration(filename=fn,
                                 d=pipeline_configuration)
