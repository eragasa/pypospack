"""
this module impplements the generalized lennard jones potential
"""
import copy
from collections import OrderedDict
import numpy as np
from pypospack.potential import PairPotential
from pypospack.potential import determine_symbol_pairs
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2019"
__license__ = "Simplified BSD License"
__version__ = 20171102

def func_cutoff_mishin2004(r, rc, hc, h0):
    ind_rc = np.ones(r.size)
    ind_rc[r > rc] = 0
    
    xrc = (r-rc)/hc
    psi_c = (xrc**4)/(1+xrc**4)

    x0 = r/h0
    psi_0 = (x0**4)/(1+x0**4)
    
    psi = psi_c * psi_0 * ind_rc

    return psi

def func_pair_generalized_lj(r,b1,b2,r1,V0,delta):
    """
    Reference:
        Y. Mishin.  Acta Materialia. 52 (2004) 1451-1467
    """

    assert isinstance(r, np.ndarray) or isinstance(r, float)
    assert isinstance(b1, float)
    assert isinstance(b2, float)
    assert isinstance(r1, float)
    assert isinstance(V0, float)
    assert isinstance(delta, float)

    z = r/r1
    assert type(z) is type(r)

    phi = (V0/(b2-b1))*((b2/(z**b1))-(b1/(z**b2)))+delta

    assert type(phi) is type(r)

    #print((5*"{:+10.4e} ").format(b1,b2,r1,V0,delta))
    #if isinstance(phi,np.ndarray):
    #    if np.isfinite(phi).any():
    #        print('need to fix here')
    return phi

def func_pair_generalized_lj_w_cutoff(r, b1, b2, r1, V0, delta, rc, hc, h0):

    phi = func_pair_generalized_lj(r,b1,b2,r1,V0,delta)
    psi = func_cutoff_mishin2004(r,rc,hc, h0)

    return psi*phi

class GeneralizedLennardJonesPotential(PairPotential):
    """Implementation of the morse potential

    Args:
        symbols(list): a list of chemical symbols.
    """

    pair_potential_parameters=['b1','b2','r1','V0','delta','rc','hc','h0']
    potential_type='general_lj'
    callable_func=func_pair_generalized_lj_w_cutoff
    def __init__(self, symbols):
        potential_type = GeneralizedLennardJonesPotential.potential_type
        PairPotential.__init__(self,
                               symbols=symbols,
                               potential_type=potential_type,
                               is_charge=False)
        self.function_pair = func_pair_generalized_lj_w_cutoff

    # this method overrides the parents stub
    def _init_parameter_names(self):
        self.symbol_pairs = list(determine_symbol_pairs(self.symbols))
        self.parameter_names = []
        for s in self.symbol_pairs:
            for p in self.pair_potential_parameters:
                self.parameter_names.append(self.PYPOSPACK_PAIR_FORMAT.format(s1=s[0], s2=s[1], p=p))
        return list(self.parameter_names)

    # this method overrides the parents stub
    def _init_parameters(self):
        self.parameters = OrderedDict()
        for v in self.parameter_names:
            self.parameters[v] = None

    # this method overrides the parent stub
    def evaluate(self, r, parameters, r_cut=None):
        """evaluate the pair potential

        This method evaluates the generalized lennard jones potential
        Args:
            r(np.ndarray): A numpy array of interatomic distances which to evaluate.
            parameters(OrderedDict): A dictionary of parameters on which to evaluate the interatomic potential.
        Notes:

            >>> symbols=['Ni']
            >>> parameters=OrderedDict([
            >>>    ('NiNi_D0',0.001114),('NiNi_a',3.429506),('NiNi_r0',2.6813)
            >>>  ])
            >>> o = MorsePotential(symbols = s)
            >>> o.evaluate(r,testing_set['parameters'])

        """
        # <----------------------------check arguments are correct
        #assert isinstance(r,np.ndarray)
        assert isinstance(parameters, OrderedDict)
        assert type(r_cut) in [int, float, type(None)]

        global_rcut = r_cut
        # <----------------------------copy a local of the parameters
        for k in self.parameters:
            self.parameters[k] = parameters[k]

        # <----------------------------evaluate the parameters now
        self.potential_evaluations = OrderedDict()
        for s in self.symbol_pairs:
            params = [self.parameters['{}_{}'.format("".join(s),k)]\
                      for k in self.pair_potential_parameters]
            self.potential_evaluations["".join(s)] \
                    = GeneralizedLennardJonesPotential.callable_func(r, *params)

        return self.potential_evaluations

    # same as parent class
    def lammps_potential_section_to_string(self):
        """needs to be overridden"""
        raise NotImplementedError

    # same as parent class
    def gulp_potential_section_to_string(self):
        """needs to be overridden"""
        raise NotImplementedError

    # same as parent class
    def phonts_potential_section_to_string(self):
        """needs to be overridden"""
        raise NotImplementedError

    # same as parent class
    def write_lammps_potential_file(self):
        """needs to be overridden"""
        raise NotImplementedError

    # same as parent class
    def write_gulp_potential_section(self):
        """needs to be overridden"""
        raise NotImplementedError

    def references(self):
        reference_dict = {}
        reference_dict['Mishin2003_NiAl_eam'] = OrderedDict()
        reference_dict['Mishin2003_NiAl_eam']['short_citation'] = "Acta Materialia 52 (2004) 1451–1467"
        reference_dict['Mishin2003_NiAl_eam']['doi'] = "doi:10.1016/j.actamat.2003.11.026"
        return reference_dict

if __name__ == "__main__":
    o = MorsePotential(symbols=['Ni'])
