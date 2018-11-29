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
pipeline_configuration[0]['function_calls'][0]['args']['cols'] = ['qoi']
pipeline_configuration[0]['function_calls'][0]['args']['clusters'] = None
pipeline_configuration[0]['function_calls'][0]['args']['kwargs'] = OrderedDict()
pipeline_configuration[0]['function_calls'][0]['args']['kwargs']['standard_scaler'] = OrderedDict()
pipeline_configuration[0]['function_calls'][0]['args']['kwargs']['standard_scaler']['with_mean'] = True
pipeline_configuration[0]['function_calls'][0]['args']['kwargs']['standard_scaler']['with_std'] = True

# define second segment (PCA transformation)
pipeline_configuration[1] = OrderedDict()
pipeline_configuration[1]['segment_type'] = 'pca'
pipeline_configuration[1]['function_calls'] = OrderedDict()
pipeline_configuration[1]['function_calls'][0]= OrderedDict()
pipeline_configuration[1]['function_calls'][0]['function'] = 'transform_pca'
pipeline_configuration[1]['function_calls'][0]['args'] = OrderedDict()
pipeline_configuration[1]['function_calls'][0]['args']['cols'] = ['n_qoi']
pipeline_configuration[1]['function_calls'][0]['args']['clusters'] = None
pipeline_configuration[1]['function_calls'][0]['args']['kwargs'] = OrderedDict()
pipeline_configuration[1]['function_calls'][0]['args']['kwargs']['pca'] = OrderedDict()
pipeline_configuration[1]['function_calls'][0]['args']['kwargs']['pca']['copy'] = False
pipeline_configuration[1]['function_calls'][0]['args']['kwargs']['pca']['whiten'] = False

# define third segment (TSNE transformation)
pipeline_configuration[2] = OrderedDict()
pipeline_configuration[2]['segment_type'] = 'manifold'
pipeline_configuration[2]['function_calls'] = OrderedDict()
pipeline_configuration[2]['function_calls'][0] = OrderedDict()
pipeline_configuration[2]['function_calls'][0]['function'] = 'transform_tsne'
pipeline_configuration[2]['function_calls'][0]['args'] = OrderedDict()
# throw out pca 0
pipeline_configuration[2]['function_calls'][0]['args']['abs_cols'] = ['pca_1', 'pca_2', 'pca_3', 'pca_4', 'pca_5', 'pca_6']
pipeline_configuration[2]['function_calls'][0]['args']['kwargs'] = OrderedDict()
pipeline_configuration[2]['function_calls'][0]['args']['kwargs']['tsne'] = OrderedDict()


if __name__ == "__main__":
    pipeline = PyposmatPipeline()
    fn = __file__.replace('.py', '.in')
    pipeline.write_configuration(filename=fn,
                                 d=pipeline_configuration)
