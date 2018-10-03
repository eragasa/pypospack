__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2018"
__license__ = "Simplified BSD License"
__version__ = 20180301

import copy
import numpy as np
from collections import OrderedDict
from pypospack.potential import PairPotential
from pypospack.potential import determine_symbol_pairs

class BornMayerPotential(PairPotential):
    """ Implementation of a Born-Mayer repulsive potential

    Args:
        symbols(list of str)
    Attributes:
        symbols(list of str)

    """
    def __init__(self,symbols):
        self.pair_potential_parameters = ['phi0','gamma','r0']
        PairPotential.__init__(self,
                symbols,
                potential_type='bornmeyer',
                is_charge=False)

    def _init_parameter_names(self):
        self.symbol_pairs = list(determine_symbol_pairs(self.symbols))
        self.parameter_names = []
        for s in self.symbol_pairs:
            for p in self.pair_potential_parameters:
                self.parameter_names.append(
                        self.PYPOSPACK_PAIR_FORMAT.format(s1=s[0],s2=s[1],p=p))
        return list(self.parameter_names)

    def _init_parameters(self):
        self.parameters = OrderedDict()
        for v in self.parameter_names:
            self.parameters[v] = None

    # this method overrides the parent stub
    def evaluate(self,r,parameters,r_cut=None):
        # <----------------------------check arguments are correct
        # assert isinstance(r,np.ndarray)
        assert isinstance(parameters,OrderedDict)
        assert type(r_cut) in [int,float,type(None)]

        # <----------------------------copy a local of the parameters 
        for k in self.parameters:
            self.parameters[k] = parameters[k]

        # <----------------------------evaluate the parameters now
        self.potential_evaluations = OrderedDict()
        for s in self.symbol_pairs:
            # <------------------------extract the parameters for symbol pair
            phi0 = self.parameters['{}{}_phi0'.format(s[0],s[1])]
            gamma  = self.parameters['{}{}_gamma'.format(s[0],s[1])]
            r0 = self.parameters['{}{}_r0'.format(s[0],s[1])]
            # <------------------------embedded morse function
            def func_bornmayer(r,phi0,gamma,r0):
                return phi0*np.exp(-gamma*(r-r0))

            _pair_name = '{}{}'.format(s[0],s[1])
            if r_cut is None:
                _V = func_bornmayer(r,phi0,gamma,r0)
                self.potential_evaluations[_pair_name] = copy.deepcopy(_V)
            else:
                _rcut = np.max(r[np.where(r < r_cut)])
                _h = r[1] - r[0]
                _V_rc = func_bornmayer(_rcut,phi0,gamma,r0)
                _V_rc_p1 = func_bornmayer(_rcut+_h,phi0,gamma,r0)
                _V_rc_m1 = func_bornmayer(_rcut-_h,phi0,gamma,r0)
                _dVdr_at_rc = (_V_rc_p1-_V_rc)/_h

                # <----- calculate morse with cutoff
                _V = func_bornmayer(r,phi0,gamma,r0)
                # <----- apply the cutoff
                _V= _V - _V_rc - _dVdr_at_rc * (r-_rcut)
                # <----- V=0, where r <= _rcut
                _V[np.where(r>=_rcut)] = 0.0
        
                self.potential_evaluations[_pair_name] = copy.deepcopy(_V)
        
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
