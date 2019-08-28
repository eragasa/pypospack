import os
import numpy as np
from numpy import linalg
from scipy import stats
import pypospack.utils
from pypospack.pyposmat.data import PyposmatConfigurationFile, PyposmatDataFile
from pypospack.statistics import GaussianKde

if __name__ == "__main__":
    data_directory = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'data','Si__sw__data','pareto_optimization_unconstrained')
    
    config_fn = os.path.join(data_directory,'pyposmat.config.in')
    data_fn = os.path.join(data_directory,'pyposmat.kde.1.out')
    
    o_config = PyposmatConfigurationFile()
    o_config.read(filename=config_fn)

    o_data = PyposmatDataFile()
    o_data.read(filename=data_fn)

    names = [k for k in o_config.qoi_names if k != "Si_dia.B"]
    X = o_data.df[o_config.qoi_names].T

    print('type(X):',type(X).__name__)
    print('X.shape:',X.shape)
    scipy_kde = stats.gaussian_kde(X)
    print('scipy_kde.d:',scipy_kde.d)
    print('scipy_kde.n:',scipy_kde.n)
    pypospack_kde = GaussianKde(X)
    print('pypospack_kde.d:',pypospack_kde.d)
    print('pypospack_kde.n:',pypospack_kde.n)

    scipy_w, scipy_v = linalg.eig(scipy_kde.covariance)
    pypospack_w, pypospack_v = linalg.eig(pypospack_kde.covariance)

    print('scipy_eig_lambda:', scipy_w)
    print('pypospack_eig_lambda:', pypospack_w)

    def eigdecomposition_cov_matrix_fix(cov):
        eig_val, eig_vec = linalg.eig(cov)
        for i,v in enumerate(eig_val):
            if v == 0:
                eig_val[i] = 1e=15
        new_cov = linalg.multi_dot(
                [eig_vec,np.diag(eig_val,linalg.inv(eig_vec))])

    pypospack_kde._data_covariance = cov
    pypospack_kde._data_inv_cov = linalg.inv(cov)


    w,v = np.linalg.eig(cov)
    print(w)
    print('----- resample points ---------------')
    n = 5
    sci_points = scipy_kde.resample(n)
    pyp_points = pypospack_kde.resample(n)

    print(sci_points.shape)
    print(pyp_points.shape)
    assert sci_points.shape == pyp_points.shape
    print(pyp_points)

    P_of_X_scipy = scipy_kde.evaluate(pyp_points)
    P_of_X_pypospack = pypospack_kde.evaluate(pyp_points)

    print(P_of_X_scipy)
    print(P_of_X_pypospack)
