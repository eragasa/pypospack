__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2017"
__license__ = "Simplified BSD License"
__version__ = 20171102

import copy
import numpy as np
from collections import OrderedDict
from pypospack.potential import PairPotential
from pypospack.potential import determine_symbol_pairs

def func_lj(r,epsilon,sigma,r_cut_pair=None):
    assert isinstance(r, np.ndarray)
    assert isinstance(epsilon,float)
    assert isinstance(sigma,float)
    assert isinstance(r_cut_pair,float) or r_cut_pair is None

    phi = 4*epsilon*((sigma/r)**12 - (sigma/r)**6)

    return phi

class LennardJonesPotential(PairPotential):
    """Potential implementation of the generalized Leonnard Jones Function

    Args:
        symbols(list): a list of chemical symbols.

    Notes:
        r_cut_global = global cutoff for Lennard Jonnes interactions [Angs]
        epsilon = controls the depth of the well
        sigma = controls the location of the well
        r_cut_pair_lj = pair cutoff for the LJ potential
        r_cut_pair_colulomb = pair cutoff for the coulomb potential
    """

    RCUT_GLOBAL_LJ_DEFAULT = 10.0
    RCUT_GLBOAL_COULOMB_DEFAULT = 10.0
    global_potential_parameters = ['r_cut_global']
    pair_potential_parameters = ['epsilon','sigma','r_cut_pair','r_cut_coulomb']

    def __init__(self,symbols):
        PairPotential.__init__(self,
                symbols,
                potential_type='lj',
                is_charge=False)


    # this method overrides the parents stub
    def _init_parameter_names(self):
        self.symbol_pairs = list(
                determine_symbol_pairs(self.symbols)
        )
        self.parameter_names = []
        for s in self.symbol_pairs:
            for p in self.global_potential_parameters:
                self.parameter_names.append(p)
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
            try:
                self.parameters[k] = parameters[k]
            except KeyError as e:
                if k == 'r_cut_global':
                    self.parameters[k] = None
                elif k.endswith('r_cut_pair'):
                    self.parameters[k] = None
                elif k.endswith('r_cut_coulomb'):
                    self.parameters[k] = None
                else:
                    raise

        # <-------------------------evaluate the parameters now
        self.potential_evaluations = OrderedDict()
        for s in self.symbol_pairs:
            _pair_name = '{}{}'.format(s[0],s[1])
            self.potential_evaluations[_pair_name] = func_lj(
                    r,
                    epsilon = self.parameters['{}{}_epsilon'.format(s[0],s[1])],
                    sigma = self.parameters['{}{}_sigma'.format(s[0],s[1])],
                    r_cut_pair = self.parameters['{}{}_r_cut_pair'.format(s[0],s[1])]
            )
        
        return self.potential_evaluations
    
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

