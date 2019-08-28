__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2017"
__license__ = "Simplified BSD License"
__version__ = 20171102

import numpy as np
from eamends_mishin2004 import func_density_mishin2004

def create_r(rmax, n):
    r = r_max/n * np.linspace(1,n,n)
    return r

if __name__ == "__main__":
    from collections import OrderedDict

    r_max = 3.
    r_N = 1000

    dens_mishin2004_parameters = OrderedDict([

