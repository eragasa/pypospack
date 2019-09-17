__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2018"
__license__ = "Simplified BSD License"
__version__ = 20180524

import copy,inspect
import numpy as np
from collections import OrderedDict
from pypospack.potential import EamEmbeddingEquationOfState

from pypospack.potential.eam_embedding_eos import (
        get_omega,
        get_pair_energy_at_a,
        get_density_at_a)

def func_zopemishin_eos(a,a0,B,E0,beta,latticetype='fcc'):
    """ return the energy value from th equation of state based on the lattice parameter

    Args:
        a (float)(np.ndarray): an array of lattice parameter values to be evaluation
        a0 (float): the equilibrium lattice constant value, in Angstroms
        B0 (float): the bulk modulus value, GPA
        E0 (float): the cohesive energy at equilibrium, in eV
        beta (float): adjustable parameter, dimensionless
        lattice_type (str): the lattice name of the structure
    Returns:
        (float)(np.ndarray): the energy implied by the equation of state evaluated at the points
        implied by the argument a

    """
    x = a/a0-1 
    Omega0 = get_omega(a0)
    alpha = np.sqrt(-(9*Omega0*B)/E0)
    E = E0 * (1+alpha*x+beta*(alpha**3)*(x**3)*(2*x+3)/((x+1)**2)) * np.exp(-alpha*x)
    return E

def func_zopemishin_embedding_function(rho,a0,B,E0,beta,latticetype='fcc'):
    """ fits the embedding function to the zope mishin equation of state

    This function only exists as a implementation prototype for optimization routines
    which requires an encapsulated function.  The EosFitter is a more general implementation,
    which should be used in most cases.

    Args:
        rho (numpy.ndarray): a list of densities for the embedding function to be evaluated at
        a0 (float): the equilibrium lattice parameters
        B0 (float): the equilibrium bulk modulus.
        beta (float): a shape parameter
        lattice_type (str): the type of lattice parameter, currently only 'fcc' is implemented

    Returns:
        (numpy.ndarray): evaluations of the embedding function for rho

    """

    func_density = o.formalisms['density']['Ni']
    func_density_param = o.parameters['popt']['density']['Ni']

    func_pair = o.formalisms['pair']['NiNi']
    func_pair_param = o.parameters['popt']['pair']['NiNi']

    func_eos = o.formalisms['embedding']['Ni']

    a_min=0
    a_max=10000
    a_tol=1.e-8

    if isinstance(rho,np.ndarray):
        rho_ = rho.tolist()
    else:
        rho_ = rho

    # calculate astars
    a = np.zeros(len(rho_))
    for i,rhostar in enumerate(rho_):
        a[i] = brentq(f=get_density_at_a,a=a_min,b=a_max,
                args=(rhostar,func_density,func_density_param),
                xtol=a_tol)
    
    E_eos = func_zopemishin_eos(a,a0,B0,E0,beta)
    E_pair = get_pair_energy_at_a(a,func_pair,func_pair_param)
    E_embedding = E_eos - E_pair

    return E_embedding

class ZopeMishinEosEmbeddingFunction(EamEmbeddingEquationOfState):
    """Implementation of the Eam Embedding class for the universal equation of state
    Args:
        symbols(list of str)
    Attributes:
        symbols(list of str)
        potential_type(str): This is set to 'eamembed_universal'
        parameter_names(list of str)
        parameters(OrderedDict): The key is the symbol associated with the
            embedding function.  On initialization, the value of each parameter
            is set to None.
        embedding(OrderedDict): The key is the symbol associated with the 
            embedding function.
        N_rho(int)
        d_rho(float)
        rho_max(float)
        rho(numpy.ndarray)
    """    
    potential_type = 'eam_embed_zopemishin'
    eos_parameter_names = ['a0', 'B', 'E0', 'beta', 'latticetype']    
    func_eos = func_zopemishin_eos
    def __init__(self,
            symbols,
            obj_density_function=None,
            obj_pair_function=None,
            lattice_type = None,
            lattice_a0 = None,
            parameters=None):

        potential_type = ZopeMishinEosEmbeddingFunction.potential_type
        EamEmbeddingEquationOfState.__init__(self,
                                             symbols=symbols,
                                             potential_type=potential_type)
        self.func_eos = ZopeMishinEosEmbeddingFunction.func_eos 
        # define member variables, and initialize with None
        self.parameters =None
        self.obj_density_fn = None
        self.obj_pair_fn = None
        self.lattice_type = None
        self.lattice_a0 = None

        # process arguments
        if parameters is not None:
            self.parameters = copy.deepcopy(parameters)

        if obj_density_function is not None:
            self.obj_density_fn = obj_density_function
        
        if obj_pair_function is not None:
            self.obj_pair_fn = obj_pair_function

        if lattice_type is not None:
            self.lattice_type = lattice_type

        if lattice_a0 is not None:
            self.lattice_a0 = lattice_a0

    def _init_parameter_names(self):
        PARAMETER_NAME_FORMAT = "{s}_{pn}"
        PARAMETER_NAMES_EA_SYMBOL = ['ecoh','latticetype','B','a0']

        self.parameter_names = []

        for s in self.symbols:
            for pn in PARAMETER_NAMES_EA_SYMBOL:
                pn_string = PARAMETER_NAME_FORMAT.format(s=s,pn=pn)
                self.parameter_names.append(pn_string)

    def _init_parameters(self):
        self.parameters = OrderedDict()
        for p in self.parameter_names:
            self.parameters[p] = None

