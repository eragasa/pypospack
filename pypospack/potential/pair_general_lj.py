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

def func_cutoff_mishin2003(r,rc,h):
    x = (r-rc)/h
    psi = np.ones(r.size) * (x**4/(1+x**4))
    
    # define the cutoff indicator, 1 except when x > x_cut
    cutoff_ind = np.ones(r.size)
    cutoff_ind[r > rc] = 0
    
    return cutoff_ind*psi

def function_generalized_lj_pair(r,b1,b2,r1,V0,delta):
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

    return phi

def func_generalized_lj_w_cutoff(r,b1,b2,r1,V0,delta,rc,h):

    phi = function_generalized_lj_pair(r,b1,b2,r1,V0,delta)
    psi = func_cutoff_mishin2003(r,rc,h)

    return psi*phi

class GeneralizedLennardJonesPotential(PairPotential):
    """Implementation of the morse potential

    Args:
        symbols(list): a list of chemical symbols.
    """


    def __init__(self, symbols):
        self.pair_potential_parameters=['b1','b2','r1','V0','delta']
        self.potential_type_name='general_lj'
        PairPotential.__init__(self,
                               symbols=symbols,
                               potential_type=self.potential_type_name,
                               is_charge=False)
        self.function_pair = function_generalized_lj_pair

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
    def evaluate(self, r, parameters, r_cut=None,h=None,mollifier_type=None):
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
            # <------------------------extract the paramters for symbol pair
            b1 = self.parameters['{}{}_b1'.format(s[0],s[1])]
            b2 = self.parameters['{}{}_b2'.format(s[0],s[1])]
            r1 = self.parameters['{}{}_r1'.format(s[0],s[1])]
            V0 = self.parameters['{}{}_V0'.format(s[0],s[1])]
            delta = self.parameters['{}{}_delta'.format(s[0],s[1])]
            parameters = (b1,b2,r1,V0,delta)
            # <------------------------determine radial cutoff
            if '{}{}_rcut'.format(s[0], s[1]) in self.parameters:
                rc = self.parameters['{}{}_rcut'.format(s[0].s[1])]
            else:
                rc = global_rcut
            # <------------------------embedded morse function
            func_pair = self.function_pair
            pair_name = '{}{}'.format(s[0], s[1])


            if r_cut is None:
                # <------------------without mollifier
                V = func_pair(r, *parameters)
                self.potential_evaluations[pair_name] = copy.deepcopy(V)
            else:
                # <-----------------derivative based mollifer
                _rcut = np.max(r[np.where(r < r_cut)])
                r_step = r[1] - r[0]

                V_rc = func_pair(_rcut, *parameters)
                V_rc_p1 = func_pair(_rcut+r_step, *parameters)
                V_rc_m1 = func_pair(_rcut-r_step, *parameters)
                dVdr_at_rc = (V_rc_p1-V_rc)/r_step

                # calculate morse with cutoff
                V = func_pair(r, *parameters)

                # apply the cutoff
                V = V - V_rc - dVdr_at_rc * (r-_rcut)

                # V=0, where r <= _rcut
                V[np.where(r >= _rcut)] = 0.0

                self.potential_evaluations[pair_name] = copy.deepcopy(V)

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
