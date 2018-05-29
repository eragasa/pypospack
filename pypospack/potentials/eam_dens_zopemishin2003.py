__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2017"
__license__ = "Simplified BSD License"
__version__ = 20180523

import copy
import numpy as np
from collections import OrderedDict
from pypospack.potential import EamDensityFunction

def func_zope_mishin_mollify(r,r_cut,h):
    """
    This is the cutoff function specified

    References:
        Zope and Mishin, Phys. Rev. B. 68, 024102 2003
    """
    x = (r-r_cut)/h
    if x >= 0.: 
        return 0.
    else:
        return x**4./(1.+x**4.)
    return psi

def func_zope_mishin_density(r,A0,B0,C0,r0,y,gamma,r_cut=None,h=None):

    _psi = None
    _rho = None

    # evaluate mollify function
    if (r_cut is not None) and (h is not None):
        _psi = func_zope_mishin_mollify(r,r_cut,h)
    else:
        _psi = 1.

    # calculate electron density.
    _rho = _psi*\
            A0*((r-r0)**y)\
            *(np.exp(-gamma*(r-r0)))\
            *(1+B0*np.exp(-gamma*(r-r0))) \
        + C0
    return _rho
    
class ZopeMishinElectronDensityFunction(EamDensityFunction):

    """
    References:
        Zope and Mishin, Phys. Rev. B. 68, 024102 2003
    """

    potential_type = 'eamdens_zopemishin2003'
    def __init__(self,symbols):
        self.density_func_parameters = [
                'A0','B0','C0','r0','y','gamma','r_cut','h'
        ]
        EamDensityFunction.__init__(
                self,
                symbols=symbols,
                potential_type='eamsdens_zopemishin2003'
        )

    # overload and implement parent method
    def _init_parameter_names(self):
        self.parameter_names = []
        for s in self.symbols:
            for p in self.density_func_parameters:
                pn = "{}_{}".format(s,p)
                self.parameter_names.append(pn)

    # overload and implement parent method
    def _init_parameters(self):
        self.parameters = OrderedDict()
        for p in self.parameter_names:
            self.parameters[p] = None

    def evaluate(self,r,parameters,r_cut=None):

        for s in self.symbols:
            for p in self.density_func_parameters:
                pn = "{}_{}".format(s,p)
                self.parameters[pn] = parameters[pn]

        for pn,pv in self.parameters.items():
            if pv is None:
                return False

        self.density_evalations = OrderedDict()
        for s in self.symbols:

            self.density_evaluations[s] = func_zope_mishin_density(
                    r = r,
                    A0 = self.parameters['{}_A0'.format(s)],
                    B0 = self.parameters['{}_B0'.format(s)],
                    C0 = self.parameters['{}_C0'.format(s)],
                    r0 = self.parameters['{}_r0'.format(s)],
                    y = self.parameters['{}_y'.format(s)],
                    gamma = self.parameters['{}_gamma'.format(s)],
                    r_cut = self.parameters['{}_r_cut'.format(s)],
                    h = self.paramters['{}_h'.format(s)]
            )
        
            retirn copy.deepcopy(self.density_evaluations)

if __name__ == '__main__':

    symbols = ['Ni']
    p = ZopeMishinElectronDensityFunction(symbols=symbols)
    print(p.potential_type)
    
