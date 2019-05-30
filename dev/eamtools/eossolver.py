__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2019"
__license__ = "MIT License"
__version__ = 20190529
# version history:
# 2019-05-29: original release of code (EJR)

import inspect
from collections import OrderedDict

import numpy as np

def determine_pair_names(symbols):
    pair_names = []
    for i1,s1 in enumerate(symbols):
        for i2,s2 in enumerate(symbols):
            if i1>=i2:
                pair_names.append("".join([s1,s2]))
    return pair_names

def get_density_at_a(a,
        rho_bar,
        func_density,
        func_density_param,
        lattice_type='fcc'):

    if lattice_type == 'fcc':
        n_NN = [12,6,24,12,24,8]
        d_NN = [a/np.sqrt(2.),a,a*np.sqrt(1.5),a*np.sqrt(2.0),a*np.sqrt(2.5),a*np.sqrt(3.0)]
    else:
        m = '{} is not an implmeented lattice_type'.format(lattice_type = 'fcc')
        raise ValueError(m)

    #calculate total density
    if isinstance(a,np.ndarray):
        rho = np.zeros(len(r))
    else:
        rho = 0.
    arg_names = [k for k in inspect.getargspec(func_density)[0] if k != 'r']
    args = [func_density_param[k] for k in arg_names]
    for k in zip(n_NN,d_NN):
        n = k[0]
        r = k[1]
        rho += n * func_density(r,*args)

    return rho - rho_bar

def get_pair_energy_at_a(a,
        func_pair,
        func_pair_param,
        lattice_type='fcc'):


    if lattice_type == 'fcc':
        n_NN = [12,6,24,12,24,8]
        d_NN = [a/np.sqrt(2.),a,a*np.sqrt(1.5),a*np.sqrt(2.0),a*np.sqrt(2.5),a*np.sqrt(3.0)]
    else:
        m = '{} is not an implmeented lattice_type'.format(lattice_type = 'fcc')

    # calculate total energy
    if isinstance(a,np.ndarray):
        phi = np.zeros(len(a))
    else:
        phi = 0.
    arg_names =[k for k in inspect.getargspec(func_pair)[0] if k != 'r']
    args = [func_pair_param[k] for k in arg_names]
    for k in zip(n_NN,d_NN):
        n = k[0]
        r = k[1]
        phi += n * func_pair(r,*args)
 
    # divide by 2, to eliminate double counting
    return 0.5 * phi

class EquationOfStateSolver(object):
    """ determines the embedding function by solving the equation of state

    """

    potential_types = ['pair','density','embedding']
    def __init__(self,symbols,r,rho):

        self.symbols = symbols
        self.r = r
        self.rho = rho
        
        self.setup_potential_dictionary(symbols=symbols)

    @property
    def pair_names(self):
        return determine_pair_names(self.symbols)

    def setup_potential_dictionary(self,symbols=None):
        if symbols is not None: self.symbols = symbols

        # set up the potential dictionary
        self.potentials = OrderedDict(
            [(k,None) for k in self.potential_types]
            )
        
        self.potentials['pair'] = OrderedDict(
            [(k,OrderedDict()) for k in self.pair_names]
            )
        for pn in self.pair_names:
            self.potentials['pair'][pn]['formalism'] = None
            self.potentials['pair'][pn]['parameters'] = None
            self.potentials['pair'][pn]['evaluation'] = None

        self.potentials['density'] = OrderedDict(
            [(k,OrderedDict()) for k in self.symbols]
        )
        for s in self.symbols:
            self.potentials['density'][s]['formalism'] = None
            self.potentials['density'][s]['parameters'] = None
            self.potentials['density'][s]['evaluation'] = None

        self.potentials['embedding'] = OrderedDict(
            [(k,OrderedDict()) for k in self.symbols]
        )
        for s in self.symbols:
            self.potentials['embedding'][s]['formalism'] = None
            self.potentials['embedding'][s]['parameters'] = None
            self.potentials['embedding'][s]['evaluation'] = None

        return self.potentials

    def set_formalism(self,func_type,name,func):
        self.potential[func_type][name] = func

    def determine_potential_type_from_token(self,token):
        if token == 'p':
            potential_type = 'pair'
        elif token == 'd':
            potential_type = 'density'
        elif token == 'e':
            potential_type = 'embedding'
        else:

            m = '{} is not a valid indicator in parameter name {}'.format(pn_token[0],pn)
            raise ValueError(m)
        return potential_type

    def set_parameters_from_dict(self,parameters):
        for pn,pv in parameters.items():
            pn_tokens = pn.split('_')
            func_type = determine_potential_type_from_token(token=pn_tokens[0])
            func_name = pn_tokens[1]
            param_name = pn_tokens[2]
             

    def get_astar(self,func_density=None,func_density_param=None,lattice_type='fcc'):
        raise NotImplementedError
