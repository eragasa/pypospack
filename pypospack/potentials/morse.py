__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2017"
__license__ = "Simplified BSD License"
__version__ = 20171102

import copy
import numpy as np
from collections import OrderedDict
from pypospack.potential import PairPotential
from pypospack.potential import determine_symbol_pairs

"""This module contains classes and test functions"""
class MorsePotential(PairPotential):

    def __init__(self,symbols):
        self.pair_potential_parameters = ['D0','a','r0']
        PairPotential.__init__(self,
                symbols,
                potential_type='morse',
                is_charge=False)

    # this method overrides the parents stub
    def _init_parameter_names(self):
        self.symbol_pairs = list(determine_symbol_pairs(self.symbols))
        self.parameter_names = []
        for s in self.symbol_pairs:
            for p in self.pair_potential_parameters:
                self.parameter_names.append(
                        self.PYPOSPACK_PAIR_FORMAT.format(s1=s[0],s2=s[1],p=p))
        return list(self.parameter_names)

    # this method overrides the parents stub
    def _init_parameters(self):
        self.parameters = OrderedDict()
        for v in self.parameter_names:
            self.parameters[v] = None

    # this method overrides the parent stub
    def evaluate(self,r,parameters,r_cut=False):
        """
        Args:
            r(numpy.ndarray): A numpy array of interatomic distances which to 
                evaluate.
            param_dict(list): A dictionary of parameters on which to evaluate
                the interatomic potential.
        """

        self.potential_evaluations = OrderedDict()

        for k in self.parameters:
            self.parameters[k] = parameters[k]

        self.potential = OrderedDict()
        for s in self.symbol_pairs:
            D0 = self.parameters['{}{}_D0'.format(s[0],s[1])]
            a  = self.parameters['{}{}_a'.format(s[0],s[1])]
            r0 = self.parameters['{}{}_r0'.format(s[0],s[1])]

            if r_cut is False:
                self.potential['{}{}'.format(s[0],s[1])] \
                        = D0 * ((1-np.exp(-a*(r-r0)))**2-1)
            else:
                idx = rcut_idx = max(np.where(r<=r_cut)[0])
                _V = D0 * ((1-np.exp(-a*(r-r0)))**2-1)
                _dr = r[1] - r[0]
                _dVdr = (_V[idx]-_V[idx-1])/_dr
                self.potential['{}{}'.format(s[0],s[1])] \
                        = _V - _dVdr*(r-r-cut)
                
        return copy.deepcopy(self.potential_evaluations)
    
    # same as parent class
    def lammps_potential_section_to_string(self):
        raise NotImplementedError
    
    # same as parent class
    def gulp_potential_section_to_string(self):
        raise NotImplementedError
    
    # same as parent class
    def phonts_potential_section_to_string(self):
        raise NotImplementedError 
    
    # same as parent class
    def write_lammps_potential_file(self):
        raise NotImplementedError

    # same as parent class
    def write_gulp_potential_section(self):
        raise NotImplementedError

    def references(self):
        reference_dict = {}
        reference_dict['LammpsMorse'] \
                = "http://lammps.sandia.gov/doc/pair_morse.html" 
        return reference_dict

if __name__ == "__main__":

    pass
