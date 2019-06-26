# -*- coding: utf-8 -*-
""" implementation of some statistical functions

This module implements some staistical tools which are not currently implemented in
any widely deployed python package.

"""
import warnings
import numpy as np
from numpy import linalg
from scipy import stats, integrate, optimize
from sklearn import preprocessing
from sklearn.exceptions import DataConversionWarning
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=DataConversionWarning)

def kde_bw_scott1992(X):
    kde = stats.gaussian_kde(X,'scott')
    return kde.factor

def kde_bw_silverman1986(X):
    kde = stats.gaussian_kde(X,'silverman')
    return kde.factor

def kde_bw_chiu1999(X):
    """
    Cross validation method of Chiu 1999

    Chiu, S.T.,  Annls. of Stat., 1991, 19, 1883-1905
    https://projecteuclid.org/download/pdf_1/euclid.aos/1176348376
    """
    def fhati(X,h,i):
        if h is float: 
            _h=h
        else: 
            _h=h[0]

        Xi = np.delete(X,i)
        kde = stats.gaussian_kde(Xi,_h)
        return kde(X[i])

    def J(X,h):
        if h is float: 
            _h = h
        else: 
            _h=h[0]

        fhat = stats.gaussian_kde(X,_h)
        #F1 = integrate.quad(lambda x: fhat(x)**2,-np.inf,np.inf)[0]
        F1 = fhat.integrate_kde(fhat)
        F2 = np.array([fhati(X,h,i) for i in range(X.shape[0])])
        return F1-2*np.mean(F2)

    #h0 = Silverman1986_h(X)
    h0 = .5
    results = optimize.minimize(lambda h: J(X,h),
                                h0,
                                method='Nelder-Mead')

    return results.x[0]

class GaussianKde(stats.gaussian_kde):

    bw_method_options = ['scott','silverman','chiu1999']

    def __init__(self, X, bw_method='scott', weights=None):
        self._initialize_scaler(X)
        X_ = self._scale_points(X)
        bw_method_ = self._get_bw_method(bw_method=bw_method)
        
        if weights is None:
            # initialize without weights
            stats.gaussian_kde.__init__(self,
                                        X_,
                                        bw_method=bw_method_)
        else:
            # initialize with weights
            stats.gaussian_kde.__init__(self,
                                        X_,
                                        bw_method=bw_method_,
                                        weights=weights)

    def _initialize_scaler(self, X):
        self.scaler = preprocessing.StandardScaler().fit(X.T)

    def _scale_points(self, X):
        X_T = self.scaler.transform(X.T)
        X_ = X_T.T

        assert X.shape == X_.shape
        return X_

    def _unscale_points(self, X):
        X_T = self.scaler.inverse_transform(X.T)
        X_ = X_T.T
        return X_

    def _get_bw_method(self,bw_method):

        bw_methods = \
            {
                'scott':'scott',
                'silverman':'silverman',
                'chiu1999':kde_bw_chiu1999
            }

        
        return bw_methods[bw_method]

    def _eigdecomposition_cov_matrix_fix(self):
        cov = self.covariance
        eig_val, eig_vec = linalg.eig(cov)
        for i,v in enumerate(eig_val):
            if v == 0:
                eig_val[i] = 1e-15 
        new_cov = linalg.multi_dot(
                [eig_vec,
                 np.diag(eig_val),
                 linalg.inv(eig_vec)])
        
        self._data_covariance = new_cov
        self._data_inv_cov = linalg.inv(new_cov)
    def evaluate(self, X):
        self._eigdecomposition_cov_matrix_fix()
        X_ = self._scale_points(X)
        return stats.gaussian_kde.evaluate(self,X_)
     
    def resample(self, size=None):
        self._eigdecomposition_cov_matrix_fix()
        X = stats.gaussian_kde.resample(self, size)
        X_ = self._unscale_points(X)
        return X_

supported_kld_distributions = [
        stats.gaussian_kde,
        GaussianKde
        ]

def kullbach_lieber_divergence(f,g,n):
    """
    Calculate the Kullbach-Lieber Divergence between f and g

    This functions does a monte carlo estimation of the Kullbach-Lieber
    divergence pre-metric between two probability distribution functions
    f and g.  

    Notes:
        for `f` and `q`, the following classes are supported. scipy.stats.kde.gaussian_kde
    Args:
        f (:obj:`scipy.stats.gaussian_kde`): A probability distribution function
        g (:obj:`scipy.stats.gaussian_kde`): A probability distribution function
        n (int): The number of sampling points

    Returns:
        tuple: returns a both the KLD convergence value, and the estimated error
        of the KLD value

    Raises:
        Run

    """
    if not any([isinstance(f,v) for v in supported_kld_distributions]):
        msg_fmt = '{} is not a supported distribution for arg f'
        msg = msg_fmt.format(type(f).__name__)
        raise TypeError(msg)
    if not any([isinstance(g,v) for v in supported_kld_distributions]):
        msg_fmt = '{} is not a supported distribution for arg g'
        msg = msg_fmt.format(type(g).__name__)
        raise TypeError(msg)
    assert isinstance(n,int)

    type_f = type(f).__name__
    type_g = type(g).__name__

    x = None # initialize, will contain x sampled from f.
    f_x = None # initialize, will contain f(x)
    g_x = None # initialize, will contain g(x)

    # draw x from f
    x = f.resample(n)

    # calculate f(x) for all x
    f_x = f.__call__(x)
    # f_x = f.evaluate(x.T)
    
    # calculate g(x) for all x
    g_x = g.__call__(x)
    # g_x = g.evaluate(x.T)
    
    with np.errstate(all='raise'):
        try:
            log_f_divide_g = np.log(f_x) - np.log(g_x)
        except FloatingPointError as e:
            for i,v in enumerate(f_x):
                if v == 0.:
                    f_x[i] = 1e-15
            for i,v in enumerate(g_x):
                if v == 0.:
                    g_x[i] = 1e-15
            log_f_divide_g = np.log(f_x) - np.log(g_x)
    # calculate the Kullbach_Lieber divergence value
    d = np.sum(f_x*log_f_divide_g)/n

    # calculate variance of the Kullbach-Lieber pre-metric
    var_d = np.var(log_f_divide_g)/n

    return d, var_d


