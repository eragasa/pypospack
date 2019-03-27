"""
this module impplements the morse potential
"""
import copy
from collections import OrderedDict
import numpy as np
from pypospack.potential import PairPotential
from pypospack.potential import determine_symbol_pairs
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2017"
__license__ = "Simplified BSD License"
__version__ = 20171102


def function_morse_potential(r, D0, a, r0):
    assert isinstance(r, np.ndarray) or isinstance(r, float)
    assert isinstance(D0, float)
    assert isinstance(a, float)
    assert isinstance(r0, float)

    return D0*(np.exp(-2*a*(r-r0))-2*np.exp(-a*(r-r0)))

class MorsePotential(PairPotential):
    """Implementation of the morse potential

    Args:
        symbols(list): a list of chemical symbols.
    """

    def __init__(self, symbols):
        self.pair_potential_parameters = ['D0', 'a', 'r0']
        PairPotential.__init__(self,
                               symbols=symbols,
                               potential_type='morse',
                               is_charge=False)

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

        This method evaluates the MorsePotential
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
            D0 = self.parameters['{}{}_D0'.format(s[0], s[1])]
            a = self.parameters['{}{}_a'.format(s[0], s[1])]
            r0 = self.parameters['{}{}_r0'.format(s[0], s[1])]

            # <------------------------determine radial cutoff
            if '{}{}_rcut'.format(s[0], s[1]) in self.parameters:
                rc = self.parameters['{}{}_rcut'.format(s[0].s[1])]
            else:
                rc = global_rcut
            # <------------------------embedded morse function
            func_morse = function_morse_potential
            #def func_morse(r,D0,a,r0):
            #    return D0*(np.exp(-2*a*(r-r0))-2*np.exp(-a*(r-r0)))

            _pair_name = '{}{}'.format(s[0], s[1])
            if r_cut is None:
                _V = func_morse(r, D0, a, r0)
                self.potential_evaluations[_pair_name] = copy.deepcopy(_V)
            else:
                _rcut = np.max(r[np.where(r < r_cut)])
                _h = r[1] - r[0]
                _V_rc = func_morse(_rcut, D0, a, r0)
                _V_rc_p1 = func_morse(_rcut+_h, D0, a, r0)
                _V_rc_m1 = func_morse(_rcut-_h, D0, a, r0)
                _dVdr_at_rc = (_V_rc_p1-_V_rc)/_h

                # <----- calculate morse with cutoff
                _V = func_morse(r, D0, a, r0)
                # <----- apply the cutoff
                _V = _V - _V_rc - _dVdr_at_rc * (r-_rcut)
                # <----- V=0, where r <= _rcut
                _V[np.where(r >= _rcut)] = 0.0

                self.potential_evaluations[_pair_name] = copy.deepcopy(_V)

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
        reference_dict['LammpsMorse'] = "http://lammps.sandia.gov/doc/pair_morse.html" 
        return reference_dict

if __name__ == "__main__":
    o = MorsePotential(symbols=['Ni'])
