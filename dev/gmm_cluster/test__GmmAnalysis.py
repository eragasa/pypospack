import pytest
import os

import numpy as np
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

max_components = 50

output_path = 'test_gmm_analysis'

def cleanup():
    if os.isdir(path):
        shutil.rmtree(output_path)

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

@pytest.mark.parametrize('names',
                        [
                            ('qois'),
                            ('parameters'),
                            ('all')
                        ])
def test__GmmAnalysis____init____arg_names_string(names):
    gmm = GmmAnalysis(configuration=o_config,
                      data=o_data,
                      names=names,
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

def table__cluster_info(gmm):
    header_row = ['cluster_id', 'N']
    for k, gmm_component in gmm.clusters.items():
        print([
            k,
            '{:5}'.format(gmm_component['N']),
            '{:10.6f}'.format(gmm_component['weight'])
        ])
        #print([k, gmm_component['N'], gmm_component['weight']])

def plot__cluster_qoi(gmm):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1,1)
    qoi_targets = gmm.configuration.qoi_targets
    n_targets = len(qoi_targets)
    for k, gmm_component in gmm.clusters.items():
        qhat = np.array(
            [gmm_component['qois']['mean'][p] for p in gmm.configuration.qoi_names]
        )
        q = np.array(
            [v for v in qoi_targets.values()]
        )
        pct_error = qhat/q - 1.
        ax.plot(range(n_targets), pct_error)
    plt.show()

def test__GmmAnalysis__do_ic_analysis():

    gmm = GmmAnalysis(configuration=o_config,
                      data=o_data,
                      names='all',
                      output_path=output_path,
                      max_components=max_components)
    gmm.do_ic_analysis()
    ic_plot_path = os.path.join(output_path,'ic_plot.png')
    assert os.path.isfile(ic_plot_path)

def table__cluster_parameters(gmm):
    header_row = ['cluster_id', ''] + [p for p in gmm.configuration.parameter_names]
    print(header_row)
    for k, gmm_component in gmm.clusters.items():
        n_parameters = len(gmm.configuration.parameter_names)
        line_format = '{:5} {:5} ' + n_parameters *'{:+5.2e}  '
        mean_info = [k, 'mu'] + [gmm_component['parameters']['mean'][p] for p in gmm.configuration.parameter_names]
        std_info = ['', 'sigma'] + [gmm_component['parameters']['std'][p] for p in gmm.configuration.parameter_names]

        print(line_format.format(*mean_info))
        print(line_format.format(*std_info))

def table__cluster_qois(gmm):

    header_row = ['cluster_id', ''] + [p for p in gmm.configuration.qoi_names]
    print(header_row)
    for k, gmm_component in gmm.clusters.items():
        n_qois = len(gmm.configuration.qoi_names)     
        line_format = '{:5} {:5} ' + n_qois *'{:+5.2e}  '
        mean_info = [k, 'mu'] + [gmm_component['qois']['mean'][p] for p in gmm.configuration.qoi_names]
        std_info = ['', 'sigma'] + [gmm_component['qois']['std'][p] for p in gmm.configuration.qoi_names]

        print(line_format.format(*mean_info))
        print(line_format.format(*std_info))

def dev__GmmAnalysis():
    max_components = 20
    n_components = 10
    gmm = GmmAnalysis(configuration=o_config,
                      data=o_data,
                      names='all',
                      output_path=output_path,
                      max_components=max_components)
    print(gmm.names)
    gmm.make_gmm_models()
    gmm.do_aic_analysis()
    gmm.do_bic_analysis()
    gmm.do_ic_analysis()
    gmm.do_cluster_analysis(n_components=n_components)

    # table__cluster_info(gmm)
    # table__cluster_parameters(gmm)
    # table__cluster_qois(gmm)
    plot__cluster_qoi(gmm)
    gmm.plot_gmm_analysis(n_components=20)

if __name__ == "__main__":
    n_components = 10
    gmm = GmmAnalysis(configuration=o_config,
                      data=o_data,
                      names='all',
                      output_path=output_path,
                      max_components=max_components)
    gmm.do_cluster_analysis(n_components=n_components)
