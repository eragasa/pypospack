import pytest
import os

import numpy as np
from sklearn.mixture import GaussianMixture

import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from manifold_analysis import ManifoldAnalysis

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

n_components = 2

output_path = 'test_pca_analysis'


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

def dev__ManifoldAnalysis():
    n_components = 2
    o = ManifoldAnalysis(configuration=o_config,
                      data=o_data,
                      names='all',
                      output_path=output_path,
                      n_components=n_components)
    # o.plot_manifold_analysis()
    o.plot_cluster_analysis(cluster_type='optics')
if __name__ == "__main__":
   dev__ManifoldAnalysis()
