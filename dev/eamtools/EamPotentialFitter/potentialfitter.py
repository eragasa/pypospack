__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2019"
__license__ = "Simplified BSD License"
__version__ = 20190624

import inspect
from collections import OrderedDict

import warnings
import numpy as np

from scipy.optimize import curve_fit

from pypospack.eamtools import SeatonSetflReader
from pypospack.eamtools import create_r
from pypospack.eamtools import create_rho
# functions for the mishin2004 density function
from pypospack.potential.eamdens_mishin2004 import (func_cutoff_mishin2004,
                                                    func_density_mishin2004,
                                                    func_density_mishin2004_w_cutoff)
# functions for the mishin2004 pair potential
from pypospack.potential.pair_general_lj import (func_pair_generalized_lj,
                                                 func_pair_generalized_lj_w_cutoff)

class EamSetflFitterException(BaseException): pass

def determine_eam_setfl_pairs(symbols):
    """ determines the order of eam pairs

    the order for eam pairs is different than for the rest of pypospack

    Args:
        symbols (list): a list of symbols

    Returns:
        (list): a list of tuples
    """

    pairs = []
    for i1, s1 in enumerate(symbols):
        for i2, s2 in enumerate(symbols):
            if i1 >= i2:
                pairs.append((s1,s2))

    return pairs

class EamPotentialFitter(object):
    """ class of determining the parameters from a SETFL file

    This class will determine the optimum parameters for a given formalism given
    a setfl file.  See the testing directory for an example on how to use this.
    Currently only works for elemental systems

    Args:
        symbols (list): a list of chemical symbols

    Attributes:
        pair_potentials (OrderedDict): key is the pair in `NiNi` format.
        density_functions (OrderedDict): key is the density in 'Ni' format.  One
            density function per element.
        embedding_functions (OrderedDict): key is the pair in `Ni' format.  One
            embedding function per element
        setfl_reader (SeatonSetflReader):
    """

    def __init__(self,symbols):
        self.symbols = symbols
        self.pair_potentials = OrderedDict([
            ("".join(k),None) for k in determine_eam_setfl_pairs(symbols)
            ])
        self.density_functions = OrderedDict([
            (k,None) for k in symbols])
        self.embedding_functions = OrderedDict([
            (k,None) for k in symbols])

        pair_names = ["".join(k) for k in determine_eam_setfl_pairs(symbols)]
        
        self.parameters = OrderedDict()
        self.parameters['p0'] = OrderedDict()
        self.parameters['p0']['pair'] = OrderedDict([(k,None) for k in pair_names])
        self.parameters['p0']['embedding'] = OrderedDict([(k,None) for k in symbols])
        self.parameters['p0']['density'] = OrderedDict([(k,None) for k in symbols])

        self.parameters['popt'] = OrderedDict()
        self.parameters['popt']['pair'] = OrderedDict([(k,None) for k in pair_names])
        self.parameters['popt']['embedding'] = OrderedDict([(k,None) for k in symbols])
        self.parameters['popt']['density'] = OrderedDict([(k,None) for k in symbols])

        self.formalisms = OrderedDict()
        self.formalisms['pair'] = OrderedDict([(k,None) for k in pair_names])
        self.formalisms['embedding'] = OrderedDict([(k,None) for k in symbols])
        self.formalisms['density'] = OrderedDict([(k,None) for k in symbols])

        self.setfl_reader = None

    def initialize_from_setfl_file(filename):
        pass
    def read_setfl_file(self,filename):
        """ reads a setfl file

        This method creates an instance of SeatonSetflReader
        Args:
            filename (str): the path to the SETFL file
        
        """
        self.setfl_reader = SeatonSetflReader(path=filename)
        self.setfl_reader.read()

    def fit_potential_pair(self,func_pair_potential,symbol_pair,param0,rlow=None):

        if not isinstance(self.setfl_reader,SeatonSetflReader):
            m = "the read_setfl_command method must be run before this method"
            raise EamSetflFitterException(m)

        pair_name = "".join(symbol_pair)

        arg_names = [k for k in inspect.signature(func_pair_potential).parameters if k!= 'r']
        p0 = [param0[k] for k in arg_names if k != 'r']

        # maximum interatomic spacing distance
        rmax = self.setfl_reader.n_r * self.setfl_reader.d_r

        # interatomic spacing distance step size
        rN = self.setfl_reader.n_r

        if rlow is None:
            r = create_r(rmax,rN)
            phi = np.array(self.setfl_reader.pair_function(pair_name))
        else:
            r_ = create_r(rmax,rN)
            phi_ = np.array(self.setfl_reader.pair_function(pair_name))
            r = r_[r_>rlow]
            phi = phi_[r_>rlow]

        # iterate until convergence
        while True:
            popt,pcov = curve_fit(func_pair_potential,r,phi,method='trf',p0=p0)

            if all([np.abs(k[1]/k[0]-1) < 0.01 for k in zip(popt,p0)]):
                break
            p0=popt

        self.parameters['p0']['pair'][pair_name] = param0
        self.parameters['popt']['pair'][pair_name] = OrderedDict(
                [(k[0],k[1]) for k in zip(arg_names,popt)])
        self.formalisms['pair'][pair_name] = func_pair_potential

    def fit_density_function(self,func_density,symbol,param0,rlow=None):

        if not isinstance(self.setfl_reader,SeatonSetflReader):
            m = "the read_setfl_command method must be run before this method"
            raise EamSetflFitterException(m)

        arg_names = [k for k in inspect.getargspec(func_density)[0] if k!= 'r']
        p0 = [param0[k] for k in arg_names]
      
        # maximum interatomic spacing distance
        rmax = self.setfl_reader.n_r * self.setfl_reader.d_r

        # interatomic spacing distance step size
        rN = self.setfl_reader.n_r
        
        if rlow is None:
            r = create_r(rmax,rN)
            rho = np.array(self.setfl_reader.density_function(symbol))
        else:
            # temporary arrays
            r_ = create_r(rmax,rN)
            rho_ = np.array(self.setfl_reader.density_function(symbol))
            # final arrays
            r = r_[r_>rlow]
            rho = rho_[r_>rlow]

        # iterate until convergence
        while True:
            popt,pcov = curve_fit(func_density,r,rho,method='trf',p0=p0)

            if all([np.abs(k[1]/k[0]-1) < 0.01 for k in zip(popt,p0)]):
                break
            p0=popt

        self.parameters['p0']['density'][symbol] = param0
        self.parameters['popt']['density'][symbol] = OrderedDict(
                    [(k[0],k[1]) for k in zip(arg_names,popt)])
        self.formalisms['density'][symbol] = func_density

    def fit_embedding_function(self,func_embedding_function,symbol,param0):
        raise NotImplementedError

    def fit_analytical_embedding_function(self):
        raise NotImplementedError

    def fit_eos_embedding_function(self,func_embedding,symbol,param0,bounds=None):

        arg_names = [
                k for k in inspect.getargspec(func_embedding)[0] \
                        if k not in ['rho','lattice_type']
        ]
        p0 = [param0[k] for k in arg_names]

        # maximum interatomic spacing distance
        rhomax = self.setfl_reader.n_rho * self.setfl_reader.d_rho

        # interatomic spacing distance step size
        rhoN = self.setfl_reader.n_rho
        
        rho = create_rho(rhomax,rhoN)
        embedding = np.array(self.setfl_reader.density_function(symbol))

        print(param0)
        # iterate until convergence
        while True:
            try:
                popt,pcov = curve_fit(func_embedding,rho,embedding,method='trf',p0=p0,bounds=bounds)
            except ValueError:
                print(param0)
                print([(k[0],k[1]) for k in zip(arg_names,popt)])

            print([(k[0],k[1]) for k in zip(arg_names,popt)])

            if all([np.abs(k[1]/k[0]-1) < 0.01 for k in zip(popt,p0)]):
                break
     
            p0=popt

        self.parameters['p0']['embedding'][symbol] = param0
        self.parameters['popt']['embedding'][symbol] = OrderedDict(
                    [(k[0],k[1]) for k in zip(arg_names,popt)])
        self.formalisms['embedding'][symbol] = func_embedding

if __name__ == "__main__":
    pass

