# -*- coding: utf-8 -*-
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2017,2018,2019"
__license__ = "Simplified BSD License"
__version__ = "1.0"

"""
This module reorgnizes the potential.py module into a separate packaage for maintainability.
- EJR, 2/22/2019
"""

import os,pathlib,copy, yaml
from collections import OrderedDict
import numpy as np
from pypospack.exceptions import BadParameterException
from pypospack.eamtools import EamSetflFile

PYPOSPACK_PAIR_FORMAT = '{s1}{s2}_{p}'
def determine_symbol_pairs(symbols):
    """determine symbol pairs

    given a list of symbols gives a list of symbol pairs in the appropriate order expected within the pypospack package.

    Args:
        symbols(list of str): a list of symbols
    Returns:
        (list of list):a list of symbol pairs
        
    """

    assert type(symbols) in [str,list]
    if type(symbols) is str: 
        symbols = [symbols]

    pairs = []
    for i1,s1 in enumerate(symbols):
        for i2,s2 in enumerate(symbols):
            if i1 <= i2:
                pairs.append([s1,s2])

    return pairs

def determine_pair_parameter_names(symbols,pair_parameter_names):
    """determine pair parameter names

    gen a list of symbols and the list of parameter names for a pair potential, this function returns the list of parameters expected for the pair potential

    Args:
        symbols(list of str): a list of symbols
        pair_parameter_names(list of str): a list of parameter_names for a pair potential.
    Returns:
        (list): a list of parameter names

    """
      

    parameter_names = []
    for s in determine_symbol_pairs(symbols):
        for p in pair_parameter_names:
            parameter_names.append(PYPOSPACK_PAIR_FORMAT.format(s1=s[0],s2=s[1],p=p))

    return parameter_names
def determine_3body_triplets(symbols):
    assert type(symbols) in [str,list]
    if type(symbols) is str:
        symbols = [symbols]

    triplets = []
    for i1,s1 in enumerate(symbols):
        for i2,s2 in enumerate(symbols):
            for i3,s3 in enumerate(symbols):
                triplets = [s1,s2,s3]

    return pairs

#------------------------------------------------------------------------------
# Base Potential
#------------------------------------------------------------------------------
from pypospack.potential.potential import Potential
from pypospack.potential.pair import PairPotential
from pypospack.potential.eam_density_function import EamDensityFunction
from pypospack.potential.eam_embedding_function import EamEmbeddingFunction
from pypospack.potential.eam_embedding_eos import EamEmbeddingEquationOfState

#------------------------------------------------------------------------------
# These are pair potentials.
#------------------------------------------------------------------------------
from pypospack.potential.pair_morse import MorsePotential
from pypospack.potential.pair_buckingham import BuckinghamPotential
from pypospack.potential.pair_bornmayer import BornMayerPotential
from pypospack.potential.pair_lj import LennardJonesPotential
from pypospack.potential.pair_general_lj import GeneralizedLennardJonesPotential
pair_potential_names = [
        'buckingham',
        'morse',
        'bornmayer',
        'lj',
        'general_lj'
]

#------------------------------------------------------------------------------
# These are three body potentials, but I don't currently have a base class prototype
# for three body potentials, so they inherent from the Potential base class
#------------------------------------------------------------------------------
from pypospack.potential.threebody_tersoff import TersoffPotential
from pypospack.potential.threebody_stillingerweber import StillingerWeberPotential

threebody_potential_names = [
        "sw",
        "tersoff"
]
#------------------------------------------------------------------------------
# These are analyical EAM density functions
#------------------------------------------------------------------------------
from pypospack.potential.eamdens_exponential import ExponentialDensityFunction
from pypospack.potential.eamdens_mishin2003 import Mishin2003DensityFunction
from pypospack.potential.eamdens_mishin2004 import Mishin2004DensityFunction
eam_density_names = [
    'eam_dens_exp',
    'eam_dens_mishin2004'
]

#------------------------------------------------------------------------------
# These are analytical EAM embedding functions
#------------------------------------------------------------------------------
from pypospack.potential.eamembed_bjs import BjsEmbeddingFunction
from pypospack.potential.eamembed_universal import UniversalEmbeddingFunction
from pypospack.potential.eamembed_fs import FinnisSinclairEmbeddingFunction

eam_embedding_analytical_names = [
        'eam_embed_bjs',
        'eam_embed_universal',
        'eam_embed_fs',
]
#------------------------------------------------------------------------------
# These potentials are EAM embedding functions which are fit by solving an
# equation of state
#------------------------------------------------------------------------------
from pypospack.potential.eamembed_eos_rose import RoseEquationOfStateEmbeddingFunction
from pypospack.potential.eamembed_eos_zopemishin import ZopeMishinEosEmbeddingFunction 

eam_embedding_eos_names = [
        'eam_embed_eos_rose',
        'eam_embed_eos_zopemishin'
]

eam_embedding_names = eam_embedding_analytical_names + eam_embedding_eos_names

from pypospack.potential.eam import EamPotential

def PotentialObjectMap(potential_type='all'):
    potential_map = OrderedDict()
    potential_map['buckingham'] = OrderedDict()
    potential_map['buckingham']['module'] = 'pypospack.potential'
    potential_map['buckingham']['class'] = 'BuckinghamPotential'

    potential_map['eam'] = OrderedDict()
    potential_map['eam']['module'] = 'pypospack.potential'
    potential_map['eam']['class'] = 'EamPotential'

    potential_map['morse'] = OrderedDict()
    potential_map['morse']['module'] = 'pypospack.potential'
    potential_map['morse']['class'] = 'MorsePotential'

    potential_map['bornmayer'] = OrderedDict()
    potential_map['bornmayer']['module'] = 'pypospack.potential'
    potential_map['bornmayer']['class'] = 'BornMayerPotential'

    potential_map['general_lj'] = OrderedDict()
    potential_map['general_lj']['module'] = 'pypospack.potential'
    potential_map['general_lj']['class'] = 'GeneralizedLennardJonesPotential'

    potential_map['tersoff'] = OrderedDict()
    potential_map['tersoff']['module'] = 'pypospack.potential'
    potential_map['tersoff']['class'] = 'TersoffPotential'

    potential_map['stillingerweber'] = OrderedDict()
    potential_map['stillingerweber']['module'] = 'pypospack.potential'
    potential_map['stillingerweber']['class'] = 'StillingerWeberPotential'

    potential_map['eam_dens_exp'] = OrderedDict()
    potential_map['eam_dens_exp']['module'] = 'pypospack.potential'
    potential_map['eam_dens_exp']['class'] = 'ExponentialDensityFunction'

    potential_map['eam_dens_mishin2003'] = OrderedDict()
    potential_map['eam_dens_mishin2003']['module'] = 'pypospack.potential'
    potential_map['eam_dens_mishin2003']['class'] = 'Mishin2003DensityFunction'

    potential_map['eam_dens_mishin2004'] = OrderedDict()
    potential_map['eam_dens_mishin2004']['module'] = 'pypospack.potential'
    potential_map['eam_dens_mishin2004']['class'] = 'Mishin2004DensityFunction'
    
    potential_map['eam_embed_bjs'] = OrderedDict()
    potential_map['eam_embed_bjs']['module'] = 'pypospack.potential'
    potential_map['eam_embed_bjs']['class'] = 'BjsEmbeddingFunction'

    potential_map['eam_embed_universal'] = OrderedDict()
    potential_map['eam_embed_universal']['module'] = 'pypospack.potential'
    potential_map['eam_embed_universal']['class'] = 'UniversalEmbeddingFunction'

    potential_map['eam_embed_fs'] = OrderedDict()
    potential_map['eam_embed_fs']['module'] = 'pypospack.potential'
    potential_map['eam_embed_fs']['class'] = 'FinnisSinclairEmbeddingFunction'

    potential_map['eam_embed_eos_rose'] = OrderedDict()
    potential_map['eam_embed_eos_rose']['module'] = 'pypospack.potential'
    potential_map['eam_embed_eos_rose']['class'] = 'RoseEquationOfStateEmbeddingFunction'

    if potential_type == 'all':
        return potential_map
    else:
        module_name = potential_map[potential_type]['module']
        class_name = potential_map[potential_type]['class']
        return module_name,class_name

# uncommented on April Fools Day, 2019
if False:
    def get_potential_map():
        """ get the potential map

        Support for interatomic potentials requires a mapping from the
        potential formalism to the class supporting the function.  Additional
        potentials to be supported can either be added here, or provided
        using any module available on the PYTHONPATH.

        Returns:
            (dict): the key is the name of the potential formalism, the value is a list where the first element contains the module. second element contains the class supporting the function

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


def func_cutoff(r,rcut,h):
    x = (r - rcut)/h
    phi = (x**4)/(1+x**4)
    # get index values where r > rcut
    r_gt_rcut = np.where(r >= rcut)
    phi[r_gt_rcut] = np.zeros(len(r_gt_rcut))
    return phi
eam_density_functions = ['eam_dens_exp']
