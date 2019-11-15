import os
import pypospack.utils
from sklearn.mixture import GaussianMixture
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

if __name__ == "__main__":

    # pypospack root directory
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

    gmm_parameters = GaussianMixture(
            n_components=10,
            covariance_type='full',
            random_state=0
    ).fit(o_data.df[o_config.parameter_names])
    o_data.df['cluster_parameter_id'] = gmm_parameters.predict(o_data.df[o_config.parameter_names])
    
    gmm_qoi = GaussianMixture(
            n_components=10,
            covariance_type='full',
            random_state=0
    ).fit(o_data.df[o_config.qoi_names])
    o_data.df['cluster_qoi_id'] = gmm_qoi.predict(o_data.df[o_config.qoi_names])

    bijective_matrix = []
    for parameter_cluster_id in range(10):
        row = []
        for qoi_cluster_id in range(10):
            row.append(len(o_data.df[ (o_data.df['cluster_parameter_id']==parameter_cluster_id) & (o_data.df['cluster_qoi_id']==qoi_cluster_id)]))
        bijective_matrix.append(row)

    for row in bijective_matrix:
        print(row)


    print(80*'-')
    # now the other way
    bijective_matrix = []
    for qoi_cluster_id in range(10):
        row = []
        for parameter_cluster_id in range(10):
            row.append(len(o_data.df[ (o_data.df['cluster_parameter_id']==parameter_cluster_id) & (o_data.df['cluster_qoi_id']==qoi_cluster_id)]))
        bijective_matrix.append(row)

    for row in bijective_matrix:
        print(row)


