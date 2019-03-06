__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2017"
__license__ = "Simplified BSD License"
__version__ = 20171102

import copy
import numpy as np
from collections import OrderedDict
from pypospack.potential import PairPotential
from pypospack.potential import determine_symbol_pairs

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

def func_generalized_lj(r,b1,b2,r1,V0,delta,r_cut=None,h=None):

    _psi = None  # mollifier function
    _phi = None  # pair potential

    # evaluate mollify function
    if (r_cut is not None) and (h is not None):
        _psi = func_zope_mishin_mollify(r,r_cut,h)
    else:
        _psi = 1.


    _phi = _psi*((V0/(b2-b1))*((b2/((r/r1)**b1))-(b1/((r/r1)**b2))) + delta)
    return _rho

class GeneralizedLennardJones(PairPotential):
    """Potential implementation of the generalized Leonnard Jones Function

    Args:
        symbols(list): a list of chemical symbols.
    """

    def __init__(self,symbols):
        self.pair_potential_parameters = [
                'b1','b2','r1','V0', 'delta','r_cut','h'
        ]
        PairPotential.__init__(self,
                symbols,
                potential_type='morse',
                is_charge=False)


    # this method overrides the parents stub
    def _init_parameter_names(self):
        self.symbol_pairs = list(
                determine_symbol_pairs(self.symbols)
        )
        self.parameter_names = []
        for s in self.symbol_pairs:
            for p in self.pair_potential_parameters:
                self.parameter_names.append(
                        self.PYPOSPACK_PAIR_FORMAT.format(
                            s1=s[0],
                            s2=s[1],
                            p=p
                        )
                )
        return list(self.parameter_names)

    # this method overrides the parents stub
    def _init_parameters(self):
        self.parameters = OrderedDict()
        for v in self.parameter_names:
            self.parameters[v] = None


    # this method overrides the parent stub
    def evaluate(self,r,parameters,r_cut=None):
        """evaluate the potential for a given set of distances for a given parameter set

        This method implements the parent method.

        Args:
            r(numpy.ndarray): A numpy array of interatomic distances which to evaluate.
            parameters(OrderedDict): A dictionary of parameters on which to evaluate the interatomic potential.
        """
        # <----------------------------check arguments are correct
        assert isinstance(r,np.ndarray)
        assert isinstance(parameters,OrderedDict)
        assert type(r_cut) in [int,float,type(None)]

        # <----------------------------copy a local of the parameters 
        for k in self.parameters:
            self.parameters[k] = parameters[k]

        # <----------------------------evaluate the parameters now
        self.potential_evaluations = OrderedDict()
        for s in self.symbol_pairs:
            _pair_name = '{}{}'.format(s[0],s[1])
            _V = func_generalized_lj(
                    r,
                    b1 = self.parameters['{}{}_b1'.format(s[0],s[1])],
                    b2 = self.parameters['{}{}_b2'.format(s[0],s[1])],
                    r1 = self.parameters['{}{}_r1'.format(s[0],s[1])],
                    V0 = self.parameters['{}{}_V0'.format(s[0],s[1])],
                    delta = self.parameters['{}{}_delta'.format(s[0],s[1])],
                    r_cut = self.parameters['{}{}_r_cut'.format(s[0],s[1])],
                    h = self.parameters['{}{}_h'.format(s[0],s[1])]
            )
            
            self.potential_evaluations[_pair_name] = copy.deepcopy(_V)
        
        return copy.deepcopy(self.potential_evaluations)
    
    # same as parent class
    def lammps_potential_section_to_string(self):
        """ lammps potential section to string """
        raise NotImplementedError
    
    # same as parent class
    def gulp_potential_section_to_string(self):
        """ lamps gulp potential section to string """
        raise NotImplementedError
    
    # same as parent class
    def phonts_potential_section_to_string(self):
        """ phonts potential section to string """
        raise NotImplementedError 
    
    # same as parent class
    def write_lammps_potential_file(self):
        """ write lammps potential section to file """
        raise NotImplementedError

    # same as parent class
    def write_gulp_potential_section(self):
        """ write gulp potential section """
        raise NotImplementedError

