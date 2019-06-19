# -*- coding: utf-8 -*-
""" implementation of some statistical functions

This module implements some staistical tools which are not currently implemented in
any widely deployed python package.

"""
import numpy as np


from scipy import stats
from sklearn import preprocessing
from sklearn.exceptions import DataConversionWarning
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=DataConversionWarning)

supported_kld_distributions = [
        stats.gaussian_kde
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

    type_f = type(f) 
    type_g = type(g)

    x = None # initialize, will contain x sampled from f.
    f_x = None # initialize, will contain f(x)
    g_x = None # initialize, will contain g(x)

    # draw x from f
    x = f.resample(n)

    # calculate f(x) for all x
    f_x = f.__call__(x)

    # calculate g(x) for all x
    g_x = g.__call__(x)

    log_f_divide_g = np.log(f_x/g_x)

    # calculate the Kullbach_Lieber divergence value
    d = np.sum(f_x*log_f_divide_g)/n

    # calculate variance of the Kullbach-Lieber pre-metric
    var_d = np.var(log_f_divide_g)/n

    return d, var_d


class GaussianKde(stats.gaussian_kde):

    def __init__(self, X, bw_method=None, weights=None):
        self.scaler = preprocessing.StandardScaler().fit(X.T)
        stats.gaussian_kde.__init__(self,self.scaler.transform(X.T).T)

    def evaluate(self, X):
        X_ = self.scaler.transform(np.array(X).T)
        return stats.gaussian_kde.evaluate(self,X_)
     
    def resample(self, size=None):
        X = stats.gaussian_kde.resample(self, size)
        X = self.scaler.inverse_transform(X.T)
        return X.T
