# -*- coding: utf-8 -*-
"""This module contains classes to interact with GULP."""
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2017"
__license__ = "Simplified BSD License"
__version__ = "1.0"

import copy, yaml
import numpy as np

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
        self.elements = []
        self.potential_type = None
        self.param_info = {}

        # for eam functions
        self.eam_pair_potential = None
        self.eam_embedding_function = None
        self.eam_density_function = None

    @property
    def symbols(self):
        return list(self.elements)

    @property
    def free_parameters(self):
        free_params = []
        for p in self.parameter_names:
            if 'equals' not in self.param_info[p]:
                free_params.append(p)

        return list(free_params)

    def read(self,fname=None):
        """ read potential information from yaml file 

        Args:
            fname(str): file to yaml file from.  If no argument is passed then
                use the filename attribute.  If the filename is set, then the
                filename attribute is also set
        """

        # set the attribute if not none
        if fname is not None:
            self.filename = fname

        try:
            pot_info = yaml.load(open(self.filename))
        except:
            raise

        # process elements of the yaml file
        self.elements = list(pot_info['elements'])
        self.parameter_names = list(pot_info['parameter_names'])
        self.potential_type = pot_info['potential_type']
        if self.potential_type == 'eam':
            self.eam_pair_potential = pot_info['eam_pair_potential']
            self.eam_embedding_function = pot_info['eam_embedding_function']
            self.eam_density_function = pot_info['eam_density_function']
        self.param_info = copy.deepcopy(pot_info['param_info'])

    def write(self,fname=None):
        """ write potential information from yaml file

        Args:
            fname(str): file to potential to.  If no argument is passed then
                use the filename attribute.  If the filename is set, then the
                filename atribute is also set.
        """

        # set the filename attribute
        if fname is not None:
            self.filename = fname

        # marshall attributes into a dict
        pot_info = {}
        pot_info['elements'] = list(self.elements)
        pot_info['parameter_names'] = list(self.parameter_names)
        pot_info['potential_type'] = self.potential_type
        if self.potential_type == 'eam':
            pot_info['eam_pair_potential'] = self.eam_pair_potential
            pot_info['eam_embedding_function'] = self.eam_embedding_function
            pot_info['eam_density_function'] = self.eam_density_function
        pot_info['param_info'] = copy.deepcopy(self.param_info)

        # dump dict to yaml
        with open(fname,'w') as f:
            yaml.dump(pot_info,f,default_flow_style=False)

    def check(self):
        """ performs sanity checks to potential configuration

        does the following checks:
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

class Potential(object):
    def __init__(self,symbols):
        self.potential_type = None
        self.symbols = list(symbols)
        self.param_names = None
        self.is_charge = None
        self.param = {}
        # self.__init_param_names()
        # self.__init_param_dict()

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

class Buckingham(Potential):

    def __init__(self,symbols):
        Potential.__init__(self,symbols)
        self.potential_type = 'buckingham'
        self.param_dict = None
        self.is_charge = True
        self._init_param_names()
        self._init_param_dict()

    @property
    def parameter_names(self): return list(self.param_names)

    def copy(self):
        clone = Buckingham(self.symbols)
        return clone

    def write_potential_file(self,fname_out,param_dict,r_cut=10.0):

        str_out = self.to_string(param_dict,r_cut)
        with open(fname_out,'w') as f:
            f.write(self.to_string(param_dict,r_cut=10.0))

    def gulp_potential_section_to_string(self,param_dict,r_cut=10.0):
        str_out = 'species\n'
        for s in self.symbols:
            str_out += "{} core {}\n".format(\
                    s,
                    param_dict['chrg_{}'.format(s)])
        str_out += 'buck\n'
        for i,si in enumerate(self.symbols):
            for j,sj in enumerate(self.symbols):
                if i<=j:
                    str_out += "{} core {} core {} {} {} {} {}\n".format(\
                            si,sj,
                            param_dict['{}{}_A'.format(si,sj)],
                            param_dict['{}{}_rho'.format(si,sj)],
                            param_dict['{}{}_C'.format(si,sj)],
                            0,
                            r_cut)
        return str_out

    def to_string(self,param_dict,r_cut=10):
        if param_dict is None:
            param_dict = copy.deepcopy(self.param_dict)
        else:
            self.param_dict = copy.deepcopy(param_dict)

        # set masses
        str_out = ''
        for i,s in enumerate(self.symbols):
            str_out += "mass {} {}\n".format(i+1,self._get_mass(s))
        str_out += "\n"

        # set groups
        for i,s in enumerate(self.symbols):
            str_out += "group {} type {}\n".format(s,i+1)
        str_out += "\n"

        # set chrg
        for i,s in enumerate(self.symbols):
            charge = self.param_dict['chrg_{}'.format(s)]
            str_out += "set group {} charge {}\n".format(s,charge)
        str_out += "\n"

        str_out += 'variable R_cut equal {}\n'.format(r_cut)
        str_out += '\n'
        str_out += 'pair_style buck/coul/long ${R_cut}\n'

        # set param values
        for i,si in enumerate(self.symbols):
            for j,sj in enumerate(self.symbols):
                if i <= j:
                    try:
                        A = self.param_dict['{}{}_A'.format(si,sj)]
                        rho = self.param_dict['{}{}_rho'.format(si,sj)]
                        C = self.param_dict['{}{}_C'.format(si,sj)]
                        str_out += "pair_coeff {} {} {} {} {} {}\n".format(\
                                i+1,j+1,A,rho,C,'${R_cut}')
                    except KeyError as ke:
                        s = str(ke)
                        print('key_requested:',s)
                        print('keys:',self.param_dict.keys())
                        raise
                    except TypeError as te: 
                        s = str(te)
                        print(self.param_dict)
                        print('key_requested:',s)
                        print('keys:',self.param_dict.keys())
                        raise
        # coulumbic charge summuation         
        str_out += "\n"
        str_out += "kspace_style pppm 1.0e-5\n"
        str_out += "\n"
        str_out += "neighbor 1.0 bin\n"
        str_out += "neigh_modify every 1 delay 0 check yes\n"


        return str_out

    def _init_param_names(self):
        """ initializes the parameter names """

        pairs = []
        for s in self.symbols:
            pairs.append('{}{}'.format(s,s))
        for i,si in enumerate(self.symbols):
            for j,sj in enumerate(self.symbols):
                if i<j:
                    pairs.append('{}{}'.format(si,sj))

        self.param_names = []
        for s in self.symbols:
            self.param_names.append('chrg_{}'.format(s))
        for p in pairs:
            self.param_names.append('{}_A'.format(p))
            self.param_names.append('{}_rho'.format(p))
            self.param_names.append('{}_C'.format(p))

    def _init_param_dict(self):
        self.param = {}
        for pn in self.param_names:
            self.param[pn] = None

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

class TersoffPotential(Potential):
    def __init__(self,symbols):
        Potential.__init__(self,symbols)
        self._pot_type = 'tersoff'
        self._determine_parameter_names()
        self._fname_potential_file = 'potential.mod'

    @property
    def parameter_names(self):
        return self._param_names

    def _determine_parameter_names(self):
        # TODO: This is only written for a single element potential
        symbols = self._symbols
        for i in range(n_symbols):
            for j in range(n_symbols):
                for k in range(n_symbols):
                    el1 = symbols[i]
                    el2 = symbols[j]
                    el3 = symbols[k]
                    self._add_parameter_names(el1,el2,el3)

    def _self_add_parameter_names(self,el1,el2,el3):
        s = "{}{}{}".format(el1,el2,el3)
        tersoff_param_names = ['m','gamma','lambda3','c','d','costheta0','n','beta',
                       'lambda2','B','R','D','lambda1','A']
        for p in param_names:
            self._param_names.append("{}_{}".format(s,p))

    def write_lammps_potential_file(self):
        fname_potential_mod = 'potential.mod'
        fname_tersoff_file = 'potential.tersoff'

        for i, el1 in enumerate(elements):
            for j, el2 in enumerate(elements):
                for k, el3 in enumerate(elements):
                    s = '{}{}{}'.format(el1,el2,el3)
                    m = self._param_dict['{}_m'.format(s)]
                    str_pot += '{} {} {}'.format(el1,el2,el3)
                    str_out += ' ' + self._param_dict['{}_gamma'.format(s)]
                    str_out += ' ' + self._param_dict['{}_lambda3'.format(s)]
                    str_out += ' ' + self._param_dict['{}_c'.format(s)]
                    str_out += ' ' + self._param_dict['{}_d'.format(s)]
                    str_out += ' ' + self._param_dict['{}_costheta0'.format(s)]
                    str_out += ' ' + self._param_dict['{}_n'.format(s)]
                    str_out += ' ' + self._param_dict['{}_beta'.format(s)]
                    str_out += ' ' + self._param_dict['{}_lambda2'.format(s)]
                    str_out += ' ' + self._param_dict['{}_B'.format(s)]
                    str_out += ' ' + self._param_dict['{}_R'.format(s)]
                    str_out += ' ' + self._param_dict['{}_D'.format(s)]
                    str_out += ' ' + self._param_dict['{}_lambda1'.format(s)]
                    str_out += ' ' + self._param_dict['{}_A'.format(s)]
                    str_out += '\n'

class EamPotential(Potential):
    def __init__(self,symbols):
        Potential.__init__(self,symbols)
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



class EamPotential(Potential):
    def __init__(self,symbols):
        Potential.__init__(self,symbols)
        self._pot_type = 'eam'
        self.param_dict = None
        self._param_names = None

        self._pair_function = None
        self._density_function = None
        self._embedding_function = None

    @property
    def pair_function(self):
        return self._pair_function

    @pair_function.setter
    def pair_function(self,func):
        self._pair_function = func

    @property
    def density_function(self):
        return self._density_function

    @density_function.setter
    def density_function(self, func):
        self._density_function = func

    @property
    def embedding_function(self):
        return self._embedding_function

    @embedding_function.setter
    def embedding_function(self, func):
        self._embedding_function = func

    @property
    def parameter_names(self):
        self._determine_parameter_names()
        return self._param_names

    def _determine_parameter_names(self):
        self._param_names  = ['p.{}'.format(v) for v in self.pair_function.parameter_names]
        self._param_names += ['d.{}'.format(v) for v in self.density_function.parameter_names]
        self._param_names += ['e.{}'.format(v) for v in self.embedding_function.parameter_names]

    #def make_setfl_file(self):

class EamElectronDensityFunction(Potential):
    pass

class ExponentialDensityFunction(EamElectronDensityFunction):
    def __init__(self,symbols):
        EamElectronDensityFunction.__init__(self,symbols)
        self._pot_type = 'elec_dens_exp'
        self.param_dict = None
        self._param_names

        self._determine_parameter_names()
        self._create_param_dictionary()

    @property
    def parameter_names(self):
        self._determine_parameter_names()
        return self._param_names

    def _create_param_dictionary(self):
        self._determine_parameter_names()
        return self._param_names

    def _determine_parameter_names(self):
        self._param_names = []
        for s in self._symbols:
            self._param_names.append('{}_rho0'.format(s))
            self._param_names.append('{}_beta'.format(s))
            self._param_names.append('{}_r0'.format(s))

    def evaluate(self,r,symbol,params,rcut=0,h=1):
        err_msg = "cannot find {} parameter for {}"

        rho0 = None
        beta = None
        r0 = None

        if '{}_rho0'.format(symbol) in self._param_names:
            rho0 = params['{}_rho0'.format(symbol)]
        else:
            str_out = err_msg.format(err_msg.format('rho0',symbol))
            raise ValueError(str_out)

        if '{}_beta'.format(symbol) in self._param_names:
            beta = params['{}_beta'.format(symbol)]
        else:
            str_out = err_msg.format(err_msg.format('beta',symbol))
            raise ValueError(str_out)

        
        if '{}_r0'.format(symbol) in self._param_names:
            r0 = params['{}_r0'.format(symbol)]
        else:
            str_out = err_msg.format(err_msg.format('r0',symbol))
            raise ValueError(str_out)

        val = rho0 * np.exp(-beta*(r/r0-1))

        if rcut != 0:
            val = val * func_cutoff(r,rcut,h)

        return val
class MorsePotential(Potential):
    def __init__(self,symbols):
        Potential.__init__(self,symbols)
        self._pot_type = 'pair_morse'
        self.param_dict = None
        self._param_names = None

        self._determine_parameter_names()
        self._create_param_dictionary()

    @property
    def parameter_names(self):
        self._determine_parameter_names()
        return self._param_names

    def _create_param_dictionary(self):
        self.param_dict = {}
        for v in self._param_names:
            self.param_dict[v] = None

    def evaluate(self,r,pair,params, rcut=0, h=1.):
        """
        Args:

            rcut(float) - the cutoff function.  If set to 0, then no cutoff
        function is applied.
        """
        err_msg = "MorsePotential cannot find {} parameter for {},{} pair"

        # free_params = De, a, re
        D0 = None
        a = None
        r0 = None

        if '{}{}_D0'.format(pair[0],pair[1]) in self._param_names:
            D0 = params['{}{}_D0'.format(pair[0],pair[1])]
        elif '{}{}_D0'.format(pair[1],pair[0]) in self._param_names:
            D0 = params['{}{}_D0'.format(pair[1],pair[0])]
        else:
            str_out = err_msg.format('D0',pair[0],pair[1])
            raise ValueError(str_out)

        if '{}{}_a'.format(pair[0],pair[1]) in self._param_names:
            a = params['{}{}_a'.format(pair[0],pair[1])]
        elif '{}{}_a'.format(pair[1],pair[0]) in self._param_names:
            a = params['{}{}_a'.format(pair[1],pair[0])]
        else:
            str_out = err_msg.format('a',pair[0],pair[1])
            raise ValueError(str_out)

        if '{}{}_r0'.format(pair[0],pair[1]) in self._param_names:
            r0 = params['{}{}_r0'.format(pair[0],pair[1])]
        elif '{}{}_r0'.format(pair[1],pair[0]) in self._param_names:
            r0 = params['{}{}_r0'.format(pair[1],pair[0])]
        else:
            str_out = err_msg.format('r0',pair[0],pair[1])
            raise ValueError(str_out)

        val = D0 * ((1 - np.exp(-a*(r-r0)))**2 -1)
        #val = D0 * (1 - np.exp(-a*(r-r0)))**2

        if rcut != 0:
            print('using rcut')
            val = val * func_cutoff(r,rcut,h)

        return val

    def _determine_parameter_names(self):
        symbols = self._symbols
        n_symbols = len(symbols)
        self._param_names = []
        for i in range(n_symbols):
            for j in range(n_symbols):
                if i <=j:
                    s_i = symbols[i]
                    s_j = symbols[j]
                    self._param_names.append("{}{}_D0".format(s_i,s_j))
                    self._param_names.append("{}{}_a".format(s_i,s_j))
                    self._param_names.append("{}{}_r0".format(s_i,s_j))

class EamEmbeddingFunction(Potential):
    def __init__(self,symbols):
        Potential.__init__(self,symbols)
        self._pot_type = None
        self.param_dict = None
        self._symbols = [s for s in symbols]

class UniversalEmbeddingFunction(EamEmbeddingFunction):
    def __init__(self,symbols):
        EamEmbeddingFunction.__init__(self,symbols)
        self._pot_type = 'embed_universal'
        self._determine_parameter_names()
        self._create_param_dictionary()

    @property 
    def parameter_names(self):
        self._determine_parameter_names()
        return self._param_names

    def _create_param_dictionary(self):
        self.param_dict = {}
        for v in self._param_names:
            self.param_dict[v] = None


    def _determine_parameter_names(self):
        self._param_names = []
        symbols = self._symbols
        for s in symbols:
            self._param_names.append("{}_F0".format(s))
            self._param_names.append("{}_p".format(s))
            self._param_names.append("{}_q".format(s))
            self._param_names.append("{}_F1".format(s))

    def eval(rho, symbol, params):

        F0 = None
        p = None
        q = None
        F1 = None

        if "{}_F0".format(symbol) in self._param_names:
            F0 = params["{}_F0".format(symbol)]
        if "{}_p".format(symbol) in self._param_names:
            p = params["{}_p".format(symbol)]
        if "{}_q".format(symbol) in self._param_names:
            q = params["{}_q".format(symbol)]
        if "{}_F1".format(symbol) in self_param_names:
            F1 = params["{}_F1".format(symbol)]

        val = F0 * ( q/(q-p)*rho**p - p/(q-p)*rho**q) + F1 *rho
        return val

class BjsEmbeddingFunction(EamEmbeddingFunction):
    def __init__(self,symbols):
        EamEmbeddingFunction.__init__(self.symbols)
        self._pot_type = 'bjs'

    @property
    def parameter_names(self):
        self._determin_parameter_names()
        return self._param_names

    def _determine_parameter_names(self):
        self._param_names = None
        sybols = self._symbols
        for s in range(symbols):
            self._param_names.append("{}_F0").format(s)
            self._param_names.append("{}_gamma").format(s)
            self._param_names.append("{}_F1").format(s)

    def eval(rho, symbol, params):

        F0 = None # free paramter
        gamma = None # free parameter
        F1 = None # free parameter

        if "{}_F0".format(symbol) in self._param_names:
            F0 = params["{}_F0".format(symbol)]
        if "{}_gamma".format(symbol) in self._param_names:
            p = params["{}_gamma".format(symbol)]
        if "{}_F1".format(symbol) in self_param_names:
            F1 = params["{}_F1".format(symbol)]

        val = F0*(1-gamma*np.ln(rho))*rho**gamma + F1*gamma
        return val

