import pytest
import os

from sklearn.mixture import GaussianMixture

import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from gmm_analysis import GmmAnalysis
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

max_components = 10

output_path = 'test_gmm_analysis'

def test__GmmAnalysis____init__():
    gmm = GmmAnalysis(configuration=o_config,
                      data=o_data,
                      names=o_config.normalized_error_names,
                      output_path=output_path,
                      max_components=max_components)
    assert isinstance(gmm.configuration, PyposmatConfigurationFile)
    assert isinstance(gmm.data, PyposmatDataFile)
    assert gmm.max_components == max_components
    assert gmm.aic_criteria is None
    assert gmm.bic_criteria is None

def test__GmmAnalysis____init____w_filenames():
    gmm = GmmAnalysis(configuration=config_fn,
                      data=data_fn,
                      names=o_config.normalized_error_names,
                      output_path=output_path,
                      max_components=max_components)
    assert isinstance(gmm.configuration, PyposmatConfigurationFile)
    assert isinstance(gmm.data, PyposmatDataFile)
    assert gmm.max_components == max_components
    assert gmm.aic_criteria is None
    assert gmm.bic_criteria is None

def test__GmmAnalysis__make_gmm_models():
    gmm = GmmAnalysis(configuration=o_config,
                      data=o_data,
                      names=o_config.normalized_error_names,
                      output_path=output_path,
                      max_components=max_components)

def test__GmmAnalysis__make_gmm_models():
    gmm = GmmAnalysis(configuration=o_config,
                      data=o_data,
                      names=o_config.normalized_error_names,
                      output_path=output_path,
                      max_components=max_components)
    gmm.make_gmm_models()
    for k in gmm.models:
        assert isinstance(gmm.models[k]['obj'], GaussianMixture)
        assert isinstance(k, int)

def test__GmmAnalysis__do_aic_analysis():
    gmm = GmmAnalysis(configuration=o_config,
                      data=o_data,
                      names=o_config.normalized_error_names,
                      output_path=output_path,
                      max_components=max_components)
    gmm.make_gmm_models()
    gmm.do_aic_analysis()
    for k in gmm.models:
        assert isinstance(gmm.models[k]['aic'], float)
    assert isinstance(gmm.aic_criteria, dict)
    assert isinstance(gmm.aic_criteria['min_components'], int)
    assert isinstance(gmm.aic_criteria['min_value'], float)

def test__GmmAnalysis__do_bic_analysis():
    gmm = GmmAnalysis(configuration=o_config,
                      data=o_data,
                      names=o_config.normalized_error_names,
                      output_path=output_path,
                      max_components=max_components)
    gmm.make_gmm_models()
    gmm.do_bic_analysis()
    for k in gmm.models:
        assert isinstance(gmm.models[k]['bic'], float)
    assert isinstance(gmm.bic_criteria, dict)
    assert isinstance(gmm.bic_criteria['min_components'], int)
    assert isinstance(gmm.bic_criteria['min_value'], float)

def test__GmmAnalysis__do_cluster_analysis():
    n_components = 10
    gmm = GmmAnalysis(configuration=o_config,
                      data=o_data,
                      names=o_config.normalized_error_names,
                      output_path=output_path,
                      max_components=max_components)
    gmm.make_gmm_models()
    gmm.do_cluster_analysis(n_components=n_components)

    for k in gmm.cluster_ids:
        assert isinstance(k, int)

def dev__GmmAnalysis__do_aic_analysis():
    gmm = GmmAnalysis(configuration=o_config,
                      data=o_data,
                      names=o_config.normalized_error_names,
                      output_path=output_path,
                      max_components=max_components)
    gmm.make_gmm_models()
    gmm.do_aic_analysis()
    for k in gmm.models:
        print(gmm.models[k])

def dev__GmmAnalysis():
    n_components = 10
    gmm = GmmAnalysis(configuration=o_config,
                      data=o_data,
                      names=o_config.normalized_error_names,
                      output_path=output_path,
                      max_components=max_components)
    gmm.make_gmm_models()
    gmm.do_aic_analysis()
    gmm.do_bic_analysis()
    gmm.do_cluster_analysis(n_components=n_components)

    from sklearn.datasets import fetch_mldata
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    names_ = gmm.names
    data_ = gmm.data.df[names_]
    qoi_tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=300)
    qoi_tsne_results = qoi_tsne.fit_transform(data_)
    gmm.data.df['qoi_tsne_2d_one'] = qoi_tsne_results[:,0]
    gmm.data.df['qoi_tsne_2d_two'] = qoi_tsne_results[:,1]
    plt.figure(figsize=(16,10))
    for cluster_id in gmm.cluster_ids:
        plt.scatter(gmm.data.df['qoi_tsne_2d_one'].loc[gmm.data.df['cluster_id'] == cluster_id],
                    gmm.data.df['qoi_tsne_2d_two'].loc[gmm.data.df['cluster_id'] == cluster_id])
    plt.show()

if __name__ == "__main__":
   dev__GmmAnalysis()

