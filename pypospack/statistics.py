# -*- coding: utf-8 -*-
""" implementation of some statistical functions

This module implements some staistical tools which are not currently implemented in
any widely deployed python package.

"""
import numpy as np
import scipy.stats

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
    if isinstance(f,scipy.stats.kde.gaussian_kde):
        raise TypeError('{} is not a supported distribution for arg f'.format(type(f)))
    if isinstance(g,scipy.stats.kde.gaussian_kde):
        raise TypeError('{} is not a supported distribution for arg g'.format(type(g)))
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


