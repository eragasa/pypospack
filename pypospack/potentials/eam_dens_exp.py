__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2017"
__license__ = "Simplified BSD License"
__version__ = 20171102

import copy
import numpy as np
from collections import OrderedDict
from pypospack.potential import EamDensityFunction

class ExponentialDensityFunction(EamDensityFunction):
    """
    Args:
        symbols(list of str)
    Attributes:
        symbols(list of str)
        potential_type(str): This is set to 'eamembed_universal'
        parameter_names(list of str)
        parameters(OrderedDict): The key is the symbol associated with the
            embedding function.  On initialization, the value of each parameter
            is set to None.
        density(OrderedDict): The key is the symbol associated with the 
            embedding function.
        N_r(int)
        d_r(float)
        r_max(float)
        r(numpy.ndarray)
    """    
    def __init__(self,symbols):
        self.density_func_parameters = ['rho0','beta','r0']
        EamDensityFunction.__init__(
                self,
                symbols=symbols,
                potential_type='eamdens_exp')

        
    def _init_parameter_names(self):
        self.parameter_names = []
        for s in self.symbols:
            for p in self.density_func_parameters:
                pn = "{}_{}".format(s,p)
                self.parameter_names.append(pn)
    
    def _init_parameters(self):
        self.parameters = OrderedDict()
        for p in self.parameter_names:
            self.parameters[p] = None

    def evaluate(self,r,parameters,r_cut=None):
        """

        Given a vector of interatomic distances, r, passed in as variable
        r, and the associated parameters of the potential.  This method
        sets the density attribute.

        Args:
            r(numpy.ndarray): This should be named as rho because it
                represents the electron density being evaluated.
            parameters(OrderedDict): This is a dictionary of the parameters
                of the embedding function for each atom.  The key is a
                string containing the ISO chemical symbol of the element.
                The value should be a numeric value.
            r_cut(float): This would be the density cutoff.  However the
                embedding energy is increasing with increasing electron
                density so the a r_cut has no physical meaning.  Any 
                variable passed into r_cut will be ignored.
        """
        assert isinstance(r,np.ndarray) or isinstance(r,float)
        assert isinstance(parameters,OrderedDict)
        assert type(r_cut) in [int,float,type(None)]
        # attribute.parameters[p] <--- arg:parameters[p]
        for s in self.symbols:
            for p in self.density_func_parameters:
                pn = "{}_{}".format(s,p)
                try:
                    self.parameters[pn] = parameters[pn]
                except KeyError as e:
                    print(80*'-')
                    print("{:^80}".format("DEBUGGING INFORMATION")
                    print(80*'-')
                    print('pn:{}'.format(pn))
                    print('arg -> parameters:')
                    for k,v in parameters.items():
                        print("    {}:{}".format(k,v))
                    print('attr -> density_func_parameters')
                    for v in self.density_func_parameters:
                        print("    {}".format(v))
                    raise 
        # cannot evaluate because
        for pn,pv in self.parameters.items():
            if pv is None:
                return False

        def func_dens_exp(r, rho0, beta,r0):
            return rho0 * np.exp(-beta*(r-r0))
        
        self.density_evaluations = OrderedDict()
        for s in self.symbols:
            rho0 = self.parameters['{}_rho0'.format(s)]
            beta = self.parameters['{}_beta'.format(s)]
            r0 = self.parameters['{}_r0'.format(s)]
            
            if r_cut is None:
                _rho = func_dens_exp(r,rho0,beta,r0)
                self.density_evaluations[s] = copy.deepcopy(_rho)
            else:
                _rho = func_dens_exp(r,rho0,beta,r0)
                _rcut = np.max(r[np.where(r<r_cut)])
                _h = r[1] - r[0]
                _rho_rc = func_dens_exp(_rcut,rho0,beta,r0)
                _rho_rc_p1 = func_dens_exp(_rcut+_h,rho0,beta,r0)
                _rho_rc_m1 = func_dens_exp(_rcut-_h,rho0,beta,r0)
                _drhodr_at_rc = (_rho_rc_p1 - _rho_rc)/_h

                _rho = _rho - _rho_rc -_drhodr_at_rc * (r-_rcut)
                _rho[np.where(r>=_rcut)] = 0
                self.density_evaluations[s] = copy.deepcopy(_rho)
              
        return copy.deepcopy(self.density_evaluations)

if __name__ == '__main__':

    symbols = ['Ni']
    p = ExponentialDensityFunction(symbols=symbols)
    print(p.potential_type)
