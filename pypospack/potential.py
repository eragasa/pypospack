# -*- coding: utf-8 -*-
"""This module contains classes to interact with GULP."""
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2017"
__license__ = "Simplified BSD License"
__version__ = "1.0"

import os,pathlib,copy, yaml
from collections import OrderedDict
import numpy as np
import pypospack.potentials

def determine_symbol_pairs(symbols):
    if not isinstance(symbols,list):
        raise ValueError("symbols must be a list")

    pairs = []
    for i1,s1 in enumerate(symbols):
        for i2,s2 in enumerate(symbols):
            if i1 <= i2:
                pairs.append([s1,s2])

    return pairs

class Potential(object):
    def __init__(self,
            symbols,
            potential_type=None,
            is_charge=None):
        self.PYPOSPACK_CHRG_FORMAT = "chrg_{s}"
        self.PYPOSPACK_PAIR_FORMAT = "{s1}{s2}_{p}"

        self.potential = None
        self.symbols = list(symbols)
        self.potential_type = potential_type
        self.is_charge = is_charge

        # these attributes will be initialized by _init_parameter_names
        self.symbol_pairs = None
        self.parameter_names = None
        self._init_parameter_names()
        
        # these attributes will be initialized by _init_parameter_names
        self.parameters = None
        self._init_parameters()

        # deprecated parameters here
        self.param = {}
        self.param_names = None         # list of str
        
    def _init_parameter_names(self):
        raise NotImplementedError
    
    def _init_parameters(self):
        raise NotImplementedError

    def evaluate(self,r,parameters,r_cut=False):
        raise NotImplementedError

    def write_lammps_potential_file(self):
        raise NotImplementedError

    def lammps_potential_section_to_string(self):
        raise NotImplementedError

    def write_gulp_potential_section(self):
        raise NotImplementedError

    def gulp_potential_section_to_string(self):
        raise NotImplementedError
    
    def _get_mass(self,element):
        if element == 'Mg':
            return 24.305
        elif element == "O":
            return 15.999
        elif element == 'Si':
            return 28.086
        elif element == 'Ni':
            return 58.6934
        else:
            raise ValueError("element {} not in database".format(element))
          
    def _get_name(self,element):
        if element == "Mg":
            return 'magnesium'
        elif element == "O":
            return 'oxygen'
        elif element == 'Si':
            return 'silicon'
        elif element == 'Ni':
            return 'nickel'
        else:
            raise ValueError('element {} not in database'.format(element))

from pypospack.potentials.tersoff import TersoffPotential

#-----------------------------------------------------------------------------
class PairPotential(Potential):
    def __init__(self,symbols,potential_type,is_charge):
        Potential.__init__(self,
                symbols=symbols,
                potential_type=potential_type,
                is_charge=is_charge)

from pypospack.potentials.morse import MorsePotential
from pypospack.potentials.buckingham import BuckinghamPotential
#-----------------------------------------------------------------------------
class EamDensityFunction(Potential):
    def __init__(self,
            symbols,
            potential_type='eamdens'):

        Potential.__init__(self,
                symbols=symbols,
                potential_type=potential_type,
                is_charge=False)

        self.density = None
from pypospack.potentials.eam_dens_exp import ExponentialDensityFunction

#------------------------------------------------------------------------------
class EamEmbeddingFunction(Potential):
    def __init__(self,
            symbols,
            potential_type='eamembed'):

        Potential.__init__(self,
                symbols=symbols,
                potential_type=potential_type,
                is_charge=False)
        self.embedding = None
from pypospack.potentials.eam_embed_bjs import BjsEmbeddingFunction
from pypospack.potentials.eam_embed_universal import UniversalEmbeddingFunction

class EamPotential(Potential):
    """
    Args:
       symbols(list of str):
       func_pairpotential(str):
       func_density(str):
       func_embedding(str):
    Attributes:
       symbols(list of str):

       obj_pair(OrderedDict of PairPotential)
       obj_density(OrderedDict of EamDensityFunction)
       obj_embedding(OrderedDict of EamEmbeddingFunction)
    """
    def __init__(self,
            symbols,
            func_pair=None,
            func_density=None,
            func_embedding=None,
            filename=None):
        
        self.obj_pair = None
        self.obj_density = None
        self.obj_embedding = None

        self.pair= None
        self.density = None
        self.embedding = None

        self.symbols = symbols

        assert type(func_pair) is str
        assert type(func_density) is str
        assert type(func_embedding) is str
        self.set_obj_pair(func_pair=func_pair)
        self.set_obj_density(func_density=func_density)
        self.set_obj_embedding(func_embedding=func_embedding)
        
        Potential.__init__(self,
                symbols=symbols,
                potential_type='eam')

    def _init_parameter_names(self):
        p_params = ['p_{}'.format(p) for p in self.obj_pair.parameter_names]
        d_params = ['d_{}'.format(p) for p in self.obj_density.parameter_names]
        e_params = ['e_{}'.format(p) for p in self.obj_embedding.parameter_names]

        self.parameter_names = list(p_params + d_params + e_params)

    def _init_parameters(self):
        self.parameters = OrderedDict()
        for p in self.parameter_names:
            self.parameters[p] = None

    def evaluate(self,r,rho,parameters,rcut=None):
        for p in self.parameters:
            self.parameters[p] = parameters[p]

        self.evaluate_pair()
        self.evaluate_density()
        self.evaluate_embedding()
    
    def set_obj_pair(self,func_pair):
        if func_pair == 'morse':
            self.obj_pair = MorsePotential(symbols=self.symbols)
        else:
            msg_err = "func_pair must be a PairPotential"
            raise ValueError(msg_err)
        if not isinstance(self.obj_pair,PairPotential):
            msg_err = "func_pair must be a PairPotential"
            raise ValueError(msg_err)

    def set_obj_density(self,func_density):
        if func_density == 'eam_dens_exp':
            self.obj_density = ExponentialDensityFunction(symbols=self.symbols)
        else:
            msg_err = "func_dens must be an EamDensityfunction"
            raise ValueError(msg_err)

        #<--- ensure that the potential configured is an EamDensityFunction
        if not isinstance(self.obj_density,EamDensityFunction):
            msg_err = (
                "func_dens must be an EamDensityFunction\n"
                 "\tfunc_density={}\n"
                 "\ttype(attr.obj_density)={}\n").format(
                         func_density,
                         type(self.obj_density))
            raise ValueError(msg_err)

    def set_obj_embedding(self,func_embedding):
        if func_embedding == 'eam_embed_universal':
            self.obj_embedding = UniversalEmbeddingFunction(symbols=self.symbols)
        elif func_embedding == 'eam_embed_bjs':
            self.obj_embedding = BjsEmbeddingFunction(symbols=self.symbols)
        else:
            msg_err = "func_embedding must be a EamEmbeddingFunction"
            raise ValueError(msg_err)

        if not isinstance(self.obj_embedding,EamEmbeddingFunction):
            msg_err = "func_embedding must be a EamEmbeddingFunction"
            raise ValueError(msg_err)

    def evaluate_pair(self,r=None,parameters=None,rcut=None):
        if parameters is not None:
            for p in self.parameters:
                self.parameters[p] = parameters[p]

        if r is not None:
            if not isinstance(r,np.ndarray):
                raise ValueError("r must be a numpy.ndarray")

        p_params = [p for p in self.parameters if p.startswith('p_')]
        p_params = [s.partition('_')[2] for s in p_params]
        self.obj_pairs.evaluate(
                r=r,
                parameters=p_params,
                r_cut=r_cut)
        self.pairpotential = copy.deepcopy(self.obj_pairs.potential)
 
    def evaluate_density(self,
            rho=None,
            parameters=None,
            rcut=None):
        
        #<--- check arguments of the function
        if parameters is not None:
            if not isinstance(parameters,dict):
                raise ValueError("parameters must be a dict")
            for p in self.parameters:
                self.parameters[p] = parameters[p]
        if rho is not None:
            if not isinstance(rho,np.ndarray):
                raise ValueError("rho must be an numpy.ndarray")

        #<--- grab the parameters of the density function
        d_params = [p for p in self.parameters if p.startswith('d_')]
        d_params = [s.partition('_')[2] for s in d_params]
        self.obj_density.evaluate(
                rho=rho,
                parameters=d_params,
                rho_cut=rho_cut)

        #<--- set the internal attribute
        self.density = copy.deepcopy(self.obj_density.density)

    def evaluate_embedding(self,
            rho=None,
            parameters=None,
            rho_cut=None):

        #<--- check arguments of the function
        if parameters is not None:
            for p in self.parameters:
                self.parameters[p] = parameters[p]
        if rho is not None:
            if isinstance(r,np.ndarray):
                self.rho = np.copy(rho)
            else:
                raise ValueError("r must be a numpy.ndarray")
            self.rho = np.copy(rho)
        if rho_cut is not None:
            if type(rho_cut) in [float,int]:
                self.rho_cut = rho_cut
            else:
                raise ValueError("rho_cut must be a numeric variable")
        #<--- grab the parameters for the embedding function
        e_params = [p for p in self.parameters if p.startswith('e_')]
        e_params = [s.partition('_')[2] for s in e_params]
        self.obj_embedding.evaluate(
                rho=self.rho,
                parameters=e_params,
                rho_cut=self.rho_cut)
        
        #<--- set the internal attribute
        self.embedding = copy.deepcopy(self.obj_embedding.embedding)
    
def PotentialObjectMap(potential_type):
    potential_map = OrderedDict()
    potential_map['buckingham'] = OrderedDict()
    potential_map['buckingham']['module'] = 'pypospack.potential'
    potential_map['buckingham']['class'] = 'BuckinghamPotential'

    potential_map['eam'] = OrderedDict()
    potential_map['eam']['module'] = 'pypospack.potential'
    potential_map['eam']['class'] = 'EmbeddedAtomModel'

    potential_map['morse'] = OrderedDict()
    potential_map['morse']['module'] = 'pypospack.potential'
    potential_map['morse']['class'] = 'MorsePotential'

    potential_map['tersoff'] = OrderedDict()
    potential_map['tersoff']['module'] = 'pypospack.potential'
    potential_map['tersoff']['class'] = 'TersoffPotential'

    potential_map['eam_dens_exp'] = OrderedDict()
    potential_map['eam_dens_exp']['module'] = 'pypospack.potential'
    potential_map['eam_dens_exp']['class'] = 'ExponentialDensityFunction'

    potential_map['eam_embed_bjs'] = OrderedDict()
    potential_map['eam_embed_bjs']['module'] = 'pypospack.potential'
    potential_map['eam_embed_bjs']['class'] = 'BjsEmbeddingFunction'

    potential_map['eam_embed_universal'] = OrderedDict()
    potential_map['eam_embed_universal']['module'] = 'pypospack.potential'
    potential_map['eam_embed_universal']['class'] = 'UniversalEmbeddingFunction'

    module_name = potential_map[potential_type]['module']
    class_name = potential_map[potential_type]['class']

    return module_name,class_name

def get_potential_map():
    """ get the potential map

    Support for interatomic potentials requires a mapping from the 
    potential formalism to the class supporting the function.  Additional
    potentials to be supported can either be added here, or provided
    using any module available on the PYTHONPATH.

    Returns:
        (dict):
            (str): name of the potential formalism
            (list): first element contains the module. second element contains
                the class supporting the function

    """
    potential_map = {\
            'buckingham':['pypospack.potential','Buckingham'],
            'eam':['pypospack.potential','EmbeddedAtomModel'],
            'morse':['pypospack.potential','MorsePotential'],
            'tersoff':['pypospack.potential','Tersoff']}
    return copy.deepcopy(potential_map)

def get_supported_potentials():
    supported_potentials = list(get_potential_map().keys())
    return supported_potentials

class PotentialInformation(object):
    """ Read/Write Empirical Interatomic Potential Information

    pypospack uses yaml files to store configuration information for required
    for optimization routines.
    
    Attributes:
        filename(str): filename of the yaml file to read/write configuration
            file.  Default is pypospack.potential.yaml'
        elements(list): list of string of chemical symbols
        parameter_names(list): list of parameter names
        potential_type(str): type of potential
        param_info(dict): param info
        eam_pair_potential(str): name of the functional form for the eam
            pair potential.  Set by default to None.
        eam_embedding_function(str): name of the functional form the eam 
            embedding function.  Set by default to None.
        eam_density_function(str): name of the functional form of the eam
            electron desnsity function.  Set by default to None.
    """
    def __init__(self):
        self.filename = 'pypospack.potential.yaml'
        self.symbols = None
        self.potential_type = None
        self.parameter_definitions = None

        # for eam functions
        self.eam_pair = None
        self.eam_embedding = None
        self.eam_density = None

        self.parameter_names = None
        self._free_parameter_names = None
    
    @property
    def free_parameter_names(self):
        self._free_parameter_names = []
        for p in self.parameter_names:
            if 'equals' not in self.parameter_definitions[p]:
                self._free_parameter_names.append(p)
        return self._free_parameter_names

    def read(self,filename=None):
        """ read potential information from yaml file 

        Args:
            fname(str): file to yaml file from.  If no argument is passed then
                use the filename attribute.  If the filename is set, then the
                filename attribute is also set
        """

        # set the attribute if not none
        if filename is not None:
            self.filename = filename

        try:
            yaml_potential_info = yaml.load(open(self.filename))
        except:
            raise

        # process elements of the yaml file
        self.symbols = list(yaml_potential_info['symbols'])
        self.parameter_names = list(yaml_potential_info['parameter_names'])
        self.potential_type = yaml_potential_info['potential_type']
        if self.potential_type == 'eam':
            self.eam_pair\
                    = yaml_potential_info['eam_pair_potential']
            self.eam_embedding\
                    = yaml_potential_info['eam_embedding_function']
            self.eam_density\
                    = yaml_potential_info['eam_density_function']
        self.parameter_definitions = copy.deepcopy(
                yaml_potential_info['parameter_definitions'])

    def write(self,
            filename=None,
            potential_configuration=None):
        """ write potential information from yaml file

        Args:
            filename(str): file to potential to.  If no argument is passed then
                use the filename attribute.  If the filename is set, then the
                filename atribute is also set.
        """

        # set the filename attribute
        if filename is not None:
            self.filename = filename

        if potential_configuration is not None:
            self.potential_configuration \
                    = copy.deepcopy(potential_configuration)
        else:
            self.potential_configuration = OrderedDict()
            self.potential_configuration['potential_type'] = self.potential_type
            self.potential_configuration['symbols'] = self.symbols
            self.potential_configuration['parameter_names'] = self.parameter_names
            if self.potential_type == 'eam':
                self.potential_configuration['eam_pair_type'] \
                        = self.eam_pair_type
                self.potential_configuration['eam_density_type'] \
                        = self.eam_density_type
                self.potential_configuration['eam_embedding_type'] \
                        = self.eam_embedding_type
            self.potential_configuration['parameter_definitions'] \
                    = copy.deepcopy(self.parameter_definitions)

        # dump dict to yaml
        with open(self.filename,'w') as f:
            yaml.dump(self.potential_configuration,
                    f,
                    default_flow_style=False)

    def check(self):
        """ performs sanity checks to potential configuration

        doe the following checks:
        1.    check to see if the potential type is supported.

        Raises:
            ValueError: if the potential type is unsupported
        """
        # initialize, set to false if a test fails
        passed_all_checks = True 

        # check1: check to see if the potential type is supporte
        if self.potential_type not in  get_supported_potentials():
            passed_all_checks = False
            raise ValueError(\
                "unsupported potential type: {}".format(self.potential_type))

            return passed_all_checks

    def get_potential_dict(self):
        potential_dict = {}
        potential_dict['potential_type'] = self.potential_type
        potential_dict['elements'] = self.elements
        potential_dict['params'] = None
        return copy.deepcopy(potential_dict)

class CutoffFunction(object):
    def __init__(self):
        pass

def cutoff_shifted_force(r,v,rcut):
    """shift potential by force

    Args:
        r (numpy.array): array of distances
        v (numpy.array): array of energies
    Returns:
        numpy.array: force shifted potential
    """

    return sv

def cutoff_shifted_energy(r,v,rcut):
    """shift potential by energy

    Args:
        r (numpy.array): array of distances
        v (numpy.array): array of energies
        rc (float): cutoff distance
    Returns:
        sv - (numpy.array) energy shifted potential
    """
    return sv

class ShiftedForceCutoff(CutoffFunction):
    def eval(self, r):
        pass

def func_cutoff(r,rcut,h):
    x = (r - rcut)/h
    phi = (x**4)/(1+x**4)
    # get index values where r > rcut
    r_gt_rcut = np.where(r >= rcut)
    phi[r_gt_rcut] = np.zeros(len(r_gt_rcut))
    return phi
