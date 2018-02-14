import numpy as np
import scipy.stats

def kullbach_lieber_divergence(f,g,n):
    """
    Calculate the Kullbach-Lieber Divergence between f and g

    This functions does a monte carlo estimation of the Kullbach-Lieber
    divergence pre-metric between two probability distribution functions
    f and g.

    Arguments:

    f - is a probability distribution function
    g - is a probability distribution function
    n - is the number of sampling points

    Returns:

    d - the Kullbach-Lieber Divergence value
    var_d - the variance of the estimate of the Kullbach_Lliber Divergence

    Notes:
    for f and g, the following classes are supported.
       scipy.stats.kde.gaussian_kde
    """
    type_f = type(f) 
    type_g = type(g)

    x = None # initialize, will contain x sampled from f.
    f_x = None # initialize, will contain f(x)
    g_x = None # initialize, will contain g(x)

    # draw x from f
    if type_f == scipy.stats.kde.gaussian_kde:
        x = f.resample(n)
    else:
        raise RuntimeError("{} isn't a supported distribution".format(type_f))

    # calculate f(x) for all x
    if type_f == scipy.stats.kde.gaussian_kde:
        f_x = f.__call__(x)
    else:
        raise RuntimeError("{} isn't a supported distribution".format(type_f))

    # calculate g(x) for all x
    if type_g == scipy.stats.kde.gaussian_kde:
        g_x = f.__call__(x)
    else:
        raise RuntimeError("{} isn't a supported distribution".format(type_g))

    log_f_divide_g = np.log(f_x/g_x)

    # calculate the Kullbach_Lieber divergence value
    d = np.sum(log_f_divide_g)/n

    # calculate variance of the Kullbach-Lieber pre-metric
    var_d = np.var(log_f_divide_g/n

    return d, var_d

