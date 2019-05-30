__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2018"
__license__ = "Simplified BSD License"
__version__ = 20180301

import copy
import numpy as np
from collections import OrderedDict
from pypospack.potential import PairPotential
from pypospack.potential import determine_symbol_pairs

class BornMayerPotential(potential.PairPotential):
    """ Implementation of a Born-Mayer repulsive potential


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
                        self.PYPOSPACK_PAIR_FORMAT.format(s1=s[0],s2=[s[1],p=p))
        return list(self.parameter_names)

    def _init_parameters(self):
        self.parameters = OrderedDict()
        for v in self.parameter_names:
            self.parameters[v] = None

    # this method overrides the parent stub
    def evaluate(self,r,parameters,r_cut=None):
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
            # <------------------------extract the parameters for symbol pair
            phi0 = self.parameters['{}{}_phi0'.format(s[0],s[1])]
            gamma  = self.parameters['{}{}_a'.format(s[0],s[1])]
            r0 = self.parameters['{}{}_r0'.format(s[0].s[1])]
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


if __name__ == "__main__":
    symbols = ['Ni']
    param_names = ['NiNi_D0', 'NiNi_a', 'NiNi_r0']
    param_dict = OrderedDict()
    param_dict['NiNi_D0'] = 0.001114
    param_dict['NiNi_a'] = 3.429506
    param_dict['NiNi_r0'] = 2.6813
    morse = potential.MorsePotential(symbols=symbols)
    assert type(morse.potential_type) is str
    assert morse.potential_type == 'morse' 
    assert type(morse.symbols) is list
    assert morse.symbols == symbols
    assert type(morse.param_names) is list
    assert morse.param_names == param_names
    assert type(morse.is_charge) is bool

    if True:
        print("potential_type:{}".format(morse.potential_type))
        print("symbols:{}".format(morse.symbols))
        print("param_names:{}".format(morse.param_names))
        print("is_charge:{}".format(morse.is_charge))

    r = r_max * np.linespace(1,100,N_r)/100
    morse.evaluate(r,param_dict)

    symbols = ['Ni','Al']
    param_names = ['NiNi_D0', 'NiNi_a', 'NiNi_r0']
    morse = potential.MorsePotential(symbols=symbols)
    assert type(morse.potential_type) is str
    assert morse.potential_type == 'morse' 
    
    print("potential_type:{}".format(morse.potential_type))
    print("symbols:{}".format(morse.symbols))
    print("param_names:{}".format(morse.param_names))

if False:
    import copy
    import matplotlib.pyplot as plt

    class PairPotentialPlot(object):
        def __init__(self,r,V,fname_plot='pair_potential.png'):
            assert isinstance(r,np.ndarray)
            assert isinstance(V,np.ndarray)

            self.r = copy.deepcopy(r)
            self.V = copy.deepcopy(V)
            self.fname_plot=fname_plot

            self.v_lim = [-0.01,0.01]
            self.r_lim = [min(self.r),max(self.r)]

        def plot(self):
            self.figure, self.ax = plt.subplots(nrows=1,ncols=1)

            self.ax.plot(self.r,self.V)
            self.ax.set_xlim(self.r_lim)
            self.ax.set_ylim(self.v_lim)

            self.figure.savefig(self.fname_plot)
            plt.close(self.figure)

            pair_plot = PairPotentialPlot(r,V) 
