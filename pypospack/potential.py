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
from pypospack.eamtools import EamSetflFile

potential_string_to_class_map = OrderedDict()
potential_string_to_class_map['bornmayer']=OrderedDict([
    ('module','pypospack.potential'),
    ('class','BornMayerPotential')])

class BadParameterException(Exception):
    def __init__(self,
            code,
            parameter_name,
            parameter_value,
            parameters):
        self.code = code
        self.parameter_name = parameter_name
        self.parameter_value = parameter_value
        self.parameters = parameters

def determine_symbol_pairs(symbols):
    if not isinstance(symbols,list):
        raise ValueError("symbols must be a list")

    pairs = []
    for i1,s1 in enumerate(symbols):
        for i2,s2 in enumerate(symbols):
            if i1 <= i2:
                pairs.append([s1,s2])

    return pairs

class BasePotential(object):
    def __init__(self,
            symbols,
            potential_type=None,
            is_charge=None):
        self.PYPOSPACK_CHRG_FORMAT = "chrg_{s}"
        self.PYPOSPACK_PAIR_FORMAT = "{s1}{s2}_{p}"
        self.PYPOSPACK_3BODY_FORMAT = "{s1}{s2}{s3}_{p}"

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
        elif element == 'Al':
            return 26.982
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
        elif element == 'Al':
            return 'aluminum'
        else:
            raise ValueError('element {} not in database'.format(element))

class Potential(object):

    def __init__(self,
            symbols,
            potential_type=None,
            is_charge=None):

        # define formatting strings
        self.PYPOSPACK_CHRG_FORMAT = "chrg_{s}"
        self.PYPOSPACK_PAIR_FORMAT = "{s1}{s2}_{p}"
        self.PYPOSPACK_3BODY_FORMAT = "{s1}{s2}{s3}_{p}"

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
        elif element == 'Al':
            return 26.982
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
        elif element == 'Al':
            return 'aluminum'
        else:
            raise ValueError('element {} not in database'.format(element))

#------------------------------------------------------------------------------
# These are three body potentials, but I don't currently have a base class prototype
# for three body potentials, so they inherent from the Potential base class
#------------------------------------------------------------------------------
from pypospack.potentials.tersoff import TersoffPotential
from pypospack.potentials.stillingerweber import StillingerWeberPotential

#-----------------------------------------------------------------------------
class PairPotential(Potential):
    def __init__(self,symbols,potential_type,is_charge):
        Potential.__init__(self,
                symbols=symbols,
                potential_type=potential_type,
                is_charge=is_charge)

        self.pair_evaluations = None

# Add the imports here for pair potentials which inherent from the PairPotential
# prototype class

from pypospack.potentials.morse import MorsePotential
from pypospack.potentials.buckingham import BuckinghamPotential
from pypospack.potentials.bornmayer import BornMayerPotential
#-----------------------------------------------------------------------------
class EamDensityFunction(Potential):
    def __init__(self,
            symbols,
            potential_type='eamdens'):

        Potential.__init__(self,
                symbols=symbols,
                potential_type=potential_type,
                is_charge=False)

        self.density_evaluations = None


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

        self.embedding_evaluations = None


from pypospack.potentials.eam_embed_bjs import BjsEmbeddingFunction
from pypospack.potentials.eam_embed_universal import UniversalEmbeddingFunction
from pypospack.potentials.eam_embed_fs import FinnisSinclairEmbeddingFunction

#------------------------------------------------------------------------------
class EamEmbeddingEquationOfState(EamEmbeddingFunction):
    def __init__(self,
            symbols,
            potential_type="eam_embed_eos"):

        EamEmbeddingFunction.__init__(self,
                symbols=symbols,
                potential_type=potential_type)

        self.is_eos = True
        self.density_fn = None
        self.pair_fn = None
        self.r_cut = None


from pypospack.potentials.eam_eos_rose import RoseEquationOfStateEmbeddingFunction
#------------------------------------------------------------------------------
from pypospack.eamtools import EamSetflFile

class EamPotential(Potential):
    """
    Args:
       symbols(list of str):
       func_pairpotential(str):
       func_density(str):
       func_embedding(str):
       filename(str):
    Attributes:
       symbols(list of str): the list of symbols
       obj_pair(OrderedDict of PairPotential)
       obj_density(OrderedDict of EamDensityFunction)
       obj_embedding(OrderedDict of EamEmbeddingFunction)
       N_r(int): number of radial points, the distance between two atoms
       r_max(float): the maximum distance between two atoms
       r_cut(float): the cutoff distance between two atoms
       N_rho(int): the number of points in the electron density evaluation
       rho_max(float): the maximum electron density
    """
    def __init__(self,
            symbols,
            func_pair=None,
            func_density=None,
            func_embedding=None,
            filename=None):

        # parameter format strings
        self.PYPOSPACK_EAM_PAIR_FORMAT = "{s1}{s2}_eam_pair_{p}"
        self.PYPOSPACK_EAM_DENSITY_FORMAT = "{s1}_eam_density_{p}"
        self.PYPOSPACK_EAM_EMBEDDING_FORMAT = "{s1}_eam_embedding_{p}"
        
        # these are pypospack.potential.Potential objects
        self.obj_pair = None
        self.obj_density = None
        self.obj_embedding = None

        self.N_r = None
        self.r_max = None
        self.r_cut = None

        self.N_rho = None
        self.rho_max = None

        # these will be numpy arrays
        self.r = None
        self.rho = None
        self.pair= None
        self.density = None
        self.embedding = None

        self.symbols = symbols

        self.setfl_filename_src = filename
        self.setfl_filename_dst = "{symbols}.eam.alloy".format(
                symbols="".join(self.symbols))
        self.setfl = None

        if filename is None:

            # If the filename is not specified, then the EamPotential must be calibrated
            # by specifying the number
            if type(func_pair) is not str:
                s = "If filename is not specified, then func_pair must be a string"
                raise ValueError(s)
            if type(func_density) is not str:
                s = "If filename is not specified, then func_density must be a string"
                raise ValueError(s)
            if type(func_embedding) is not str:
                s = "If filename is not specified, then func_embeddding must be a string"
                raise ValueError(s)

            self.set_obj_pair(func_pair=func_pair)
            self.set_obj_density(func_density=func_density)
            self.set_obj_embedding(func_embedding=func_embedding)
        else:
            raise NotImplementedError()

        Potential.__init__(self,
                symbols=symbols,
                potential_type='eam',
                is_charge = False)

    def lammps_potential_section_to_string(self,setfl_dst_filename):

        # set masses
        str_out = ''
        for i,s in enumerate(self.symbols):
            str_out += "mass {} {}\n".format(i+1,self._get_mass(s))
        str_out += "\n"

        # set groups
        for i,s in enumerate(self.symbols):
            str_out += "group {} type {}\n".format(s,i+1)
        str_out += "\n"

        str_out += "pair_style eam/alloy\n"
        str_out += "pair_coeff * * {setfl_dst_filename} {str_symbols}\n".format(
                    setfl_dst_filename=setfl_dst_filename,
                    str_symbols=" ".join(self.symbols))

        return str_out


    def _init_parameter_names(self):
        if all([self.obj_pair is not None,
                self.obj_density is not None,
                self.obj_embedding is not None]):
            p_params = ['p_{}'.format(p) for p in self.obj_pair.parameter_names]
            d_params = ['d_{}'.format(p) for p in self.obj_density.parameter_names]
            e_params = ['e_{}'.format(p) for p in self.obj_embedding.parameter_names]
            self.parameter_names = list(p_params + d_params + e_params)
        else:
            self.parameter_names = None

    def _init_parameters(self):
        if self.parameter_names is None:
            return

        self.parameters = OrderedDict()
        for p in self.parameter_names:
            self.parameters[p] = None

    def determine_r_max(self,a0,latt_type):
        _a0 = a0

        if latt_type == 'fcc':
            _d_1NN = 0.707 * _a0
            _d_2NN = 1.000 * _a0
            _d_3NN = 1.225 * _a0
            _d_4NN = 1.414 * _a0
            _d_5NN = 1.581 * _a0

        # r_max should be between the 3NN and 4NN
        _rcut = 0.5 * (_d_3NN + _d_4NN)
        return _rcut

    def determine_rho_max(self,a0,latt_type):
        _a0 = a0

        #extract parameters for the density function
        _parameters = OrderedDict()
        for k,v in self.parameters.items():
            if k.startswith('d_'):
                _parameter_name = k[2:]
                _parameters[_parameter_name] = v

        if latt_type == 'fcc':
            _d_1NN = 0.707 * _a0
            _d_2NN = 1.000 * _a0
            _d_3NN = 1.225 * _a0
            _d_O = 0.866 * _a0
            _d_T = 0.433 * _a0

            _natoms_1NN = 12
            _natoms_2NN = 6
            _natoms_O = 4
            _natoms_T = 2

        s = 'Ni'
        _rhomax = _natoms_1NN * self.obj_density.evaluate(_d_1NN,_parameters)[s]
        _rhomax += _natoms_2NN * self.obj_density.evaluate(_d_2NN,_parameters)[s]
        _rhomax += _natoms_O * self.obj_density.evaluate(_d_O,_parameters)[s]
        _rhomax += _natoms_T * self.obj_density.evaluate(_d_T,_parameters)[s]

        return float(_rhomax)

    def write_setfl_file(self,filename,symbols,
            Nr,rmax,rcut,
            Nrho,rhomax,
            parameters):
        assert type(filename) is str
        assert type(Nr) is int
        assert type(rmax) in [int,float]
        assert type(rcut) in [int,float]
        assert type(Nrho) is int
        assert type(rhomax) in [int,float]

        r = rmax * np.linspace(1,Nr,Nr)/Nr
        rho = rhomax * np.linspace(1,Nrho,Nrho)/Nrho

        self.evaluate(
                r=r,
                rho=rho,
                rcut=rcut,
                parameters=parameters)

        setfl_file = EamSetflFile()
        setfl_file.write(
                filename=filename,
                symbols=symbols,
                r=self.r,
                rho=self.rho,
                rcut=self.r_cut,
                pair=self.pair,
                density=self.density,
                embedding=self.embedding)

    def evaluate(self,r,rho,rcut,parameters):
        assert isinstance(r,np.ndarray)
        assert isinstance(rho,np.ndarray)
        assert type(rcut) in [float,int,type(None)]
        assert type(parameters) in [OrderedDict,dict]

        for p in self.parameters:
            self.parameters[p] = parameters[p]

        self.r_cut = rcut
        self.r = copy.deepcopy(r)
        self.rho = copy.deepcopy(rho)
        self.evaluate_pair(r=r,rcut=rcut)
        self.evaluate_density(r=r,rcut=rcut)
        self.evaluate_embedding(rho)

    def set_obj_pair(self,func_pair):
        if func_pair == 'morse':
            self.obj_pair = MorsePotential(symbols=self.symbols)
        elif func_pair == 'bornmayer':
            self.obj_pair = BornMayerPotential(symbols=self.symbols)
        else:
            msg_err = ["func_pair must be a PairPotential"]
            msg_err = ["type(func_pair)={}".format(str(type(func_pair)))]
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

        # this programming is a little sloppy.  some software architecture needs to be implemented here
        # for a more clean object oriented approach
        if func_embedding == 'eam_embed_universal':
            self.obj_embedding = UniversalEmbeddingFunction(symbols=self.symbols)
        elif func_embedding == 'eam_embed_bjs':
            self.obj_embedding = BjsEmbeddingFunction(symbols=self.symbols)
        elif func_embedding == 'eam_embed_fs':
            self.obj_embedding = FinnisSinclairEmbeddingFunction(symbols=self.symbols)
        elif func_embedding == 'eam_embed_eos_rose':
            self.obj_embedding = RoseEquationOfStateEmbeddingFunction(symbols=self.symbols)
        else:
            msg_err = "func_embedding must be a EamEmbeddingFunction"
            raise ValueError(msg_err)

        # If the Embedding Function is determined from the Equation of State, then
        # the density function and the pair function need to be provided to the
        # embedding function object since these values are determine by numerical
        # inversion of the equation of state

        if isinstance(self.obj_embedding,EamEmbeddingEquationOfState):
            self.density_fn = self.obj_density 
            self.pair_fn = self.obj_pair 

        if not isinstance(self.obj_embedding,EamEmbeddingFunction):
            msg_err = "func_embedding must be a EamEmbeddingFunction"
            raise ValueError(msg_err)

    def evaluate_pair(self,r,parameters=None,rcut=None):
        assert isinstance(r,np.ndarray)

        if parameters is not None:
            for p in self.parameters:
                self.parameters[p] = parameters[p]

        # pair potential parameters are prepended with 'p_'
        p_params = [p for p in self.parameters if p.startswith('p_')]

        # subselect the parameters required for the pair potential
        _parameters = OrderedDict()
        for p in p_params:
            _parameter_name = p.partition('_')[2]
            _parameters[_parameter_name] = self.parameters[p]

        self.obj_pair.evaluate(
                r=r,
                parameters=_parameters,
                r_cut=rcut)

        self.pair = copy.deepcopy(self.obj_pair.potential_evaluations)

    def evaluate_density(self,
            r,
            rcut=None,
            parameters=None):
        assert isinstance(r,np.ndarray) or isinstance(r,float)
        assert isinstance(rcut,float) or rcut is None
        assert isinstance(parameters,dict) or parameters is None

        #<--- check arguments of the function
        if parameters is not None:
            for p in self.parameters:
                self.parameters[p] = parameters[p]

        #<--- grab the parameters of the density function
        d_params = [p for p in self.parameters if p.startswith('d_')]

        _parameters = OrderedDict()
        for p in d_params:
            _parameter_name = p.partition('_')[2]
            _parameters[_parameter_name] = self.parameters[p]

        _dens_eval= self.obj_density.evaluate(
                r=r,
                parameters=_parameters,
                r_cut=rcut)

        #<--- set the internal attribute
        self.density = copy.deepcopy(_dens_eval)

        return _dens_eval

    def evaluate_embedding(self,
            rho=None,
            parameters=None):
        #<--- check arguments of the function
        if parameters is not None:
            for p in self.parameters:
                self.parameters[p] = parameters[p]
        if rho is not None:
            if isinstance(rho,np.ndarray):
                self.rho = np.copy(rho)
            else:
                raise ValueError("r must be a numpy.ndarray")
            self.rho = np.copy(rho)
        #<--- grab the parameters for the embedding function
        e_params = [p for p in self.parameters if p.startswith('e_')]

        _parameters = OrderedDict()
        for p in e_params:
            _parameter_name = p.partition('_')[2]
            _parameters[_parameter_name] = self.parameters[p]

        #<--- for solving embedding function implicitly from the RoseEquationOfState
        if type(self.obj_embedding) == RoseEquationOfStateEmbeddingFunction:
            if parameters is not None:
                self.obj_embedding.evaluate(
                    rho=self.rho,
                    parameters=parameters, # pass in all parameters not just a subset
                    o_pair=self.obj_pair,
                    o_density=self.obj_density)
            else:
                self.obj_embedding.evaluate(
                    rho=self.rho,
                    r=self.r,
                    parameters=self.parameters, # pass in all parameters not just a subset
                    o_pair=self.obj_pair,
                    o_density=self.obj_density)

        else:
            self.obj_embedding.evaluate(
                rho=self.rho,
                parameters=_parameters)

        #<--- set the internal attribute
        self.embedding = copy.deepcopy(self.obj_embedding.embedding_evaluations)

def PotentialObjectMap(potential_type):
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

    potential_map['tersoff'] = OrderedDict()
    potential_map['tersoff']['module'] = 'pypospack.potential'
    potential_map['tersoff']['class'] = 'TersoffPotential'

    potential_map['stillingerweber'] = OrderedDict()
    potential_map['stillingerweber']['module'] = 'pypospack.potential'
    potential_map['stillingerweber']['class'] = 'StillingerWeberPotential'

    potential_map['eam_dens_exp'] = OrderedDict()
    potential_map['eam_dens_exp']['module'] = 'pypospack.potential'
    potential_map['eam_dens_exp']['class'] = 'ExponentialDensityFunction'

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


def func_cutoff(r,rcut,h):
    x = (r - rcut)/h
    phi = (x**4)/(1+x**4)
    # get index values where r > rcut
    r_gt_rcut = np.where(r >= rcut)
    phi[r_gt_rcut] = np.zeros(len(r_gt_rcut))
    return phi
