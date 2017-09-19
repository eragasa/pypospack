# -*- coding: utf-8 -*-
"""This module contains the pyposmat engine for parameterization"""
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2017"
__license__ = "Simplified BSD License"
__version__ = "1.0"

import time, yaml, copy
import os, shutil, subprocess

import numpy as np
import scipy.stats

def get_supported_qois():
    supported_qois = ['a0','a1','a2','a3',
                     'alpha','beta','gamma',
                     'c11','c12','c44',
                     'bulk_modulus',
                     'shear_modulus',
                     'defect_energy',
                     'surface_energy',
                     'stacking_fault_energy',
                     'total_energy']
    return supported_qois

class EipFittingError(Exception):
    pass

class EipFittingEngine(object):
    """ Generic Fitting Engine

    This fitting engine does not have an algorithm.

    Args:
        fname_config_pyposmat(string): filename of the configuration file.
           default is pyposmat.config
        fname_config_potential(string): filename of the potential file.
           default is pyposmat.potential
        fname_config_qoi(string): filename of the qoi file
        random_seed(int): random seed to use.  Default is set to None, which
           generates the random seed automatically.
        restart(bool): when set to True, attempts to restart simulations from
           existing information contained in the directory

    Attributes:
        fname_config_potential(str)
        fname_config_qoi(str)
        random_seed(int)
        restart(bool)
    """
    def __init__(self,
            fname_config_structures = "pypospack.structure.yaml",
            fname_config_potential = "pypospack.potential.yaml",
            fname_config_qoi = "pypospack.qoi.yaml",
            fname_results = "results.out",
            fname_log = "pypospack.log",
            random_seed = None,restart = False):

        self.supported_qoi = list(get_supported_qois())

        self.fname_config_structures = fname_config_structures
        self.fname_config_potential = fname_config_potential
        self.fname_config_qoi = fname_config_qoi
        
        self.restart = restart
        if self.restart is True:
            raise NotImplementedError("Restart method not implemented")

        self._set_random_seed(random_seed)    

        # determine output
        self._configure_results_file(fname_results)
        self._configure_log_file(fname_log)

        # read configuration files
        self.structures = StructureDatabase(self.fname_config_structures)
        self.potential = PotentialInformation(self.fname_config_potential)
        self.qois = QoiDatabase(self.fname_config_qoi)
    def _set_random_seed(self,seed):
        # set the random seed
        np.random.seed(seed)
        # get the random seed from numpy
        self.random_seed = np.random.get_state()[1][0]

    def _configure_results_file(self,fname):
        self._f_results = open(fname,'w')

    def _configure_log_file(self,fname):
        self._f_log = open(fname,'w')

    def _configure_potential(self):
        self._log('configure the potential')

class StructureDatabase(object):
    """ structure database 

    Attributes:
        filename(str): file to read/write the yaml file
        directory(str): the directory of the structure database
        structures(dict): key is structure name, value is another dict with 
            key/value pairs for 'filename' and 'filetype'
    """

    def __init__(self):
        self.filename = 'pypospack.structure.yaml'
        self.directory = None
        self.structures = {}

    def add_structure(self,name,filename,filetype):
        self.structures[name] = {'filename':filename,'filetype':filetype}

    def read(self,fname = None):
        """ read qoi configuration from yaml file

        Args:
            fname(str): file to read yaml file from.  If no argument is passed 
                then use the filename attribute.  If the filename is set, then 
                the filename attribute is also set.
        """

        # set the attribute if not none
        if fname is not None:
            self.filename = fname

        try:
            self.structure_db = yaml.load(open(self.filename))
        except:
            raise

        self.directory = self.structure_db['directory']
        self.structures = copy.deepcopy(self.structure_db['structures'])

    def write(self,fname = None):
        if fname is not None:
            self.filename = fname

        # marshall attributes into a dictionary
        self.structure_db = {}
        self.structure_db['directory'] = self.directory
        self.structure_db['structures'] = {}
        self.structure_db['structures'] = copy.deepcopy(self.structures)

        # dump to as yaml file
        with open(fname,'w') as f:
            yaml.dump(self.structure_db, f, default_flow_style=False)

class QoiDatabase(object):
    """ Qoi Database 
    
        Attributes:
            filename(str): file to read/write to yaml file
            qois(dict): key is qoi name, value contains a variety of values
    """
        
    def __init__(self):
        self.filename = 'pypospack.qoi.yaml'
        self.qois = {}

    def add_qoi(self,name,qoi,structures,target):
        """ add a qoi

        Args:
            name(str): name of the qoi.  Usually <structure>.<qoi>.
            qoi(str): name of the qoi.
            structures(list): list of structures
        """

        # make a copy of structures
        _structures = None
        if isinstance(structures,str):
            _structures = [structures]
        else:
            _structures = list(structures)

        self.qois[name] = {\
                'qoi':qoi,
                'structures':copy.deepcopy(structures),
                'target':target}

    def read(self,fname=None): 
        """ read qoi configuration from yaml file

        Args:
            fname(str): file to read yaml file from.  If no argument is passed
                then use the filename attribute.  If the filename is set, then
                the filename attribute is also set.
        """

        # set the attribute if not none
        if fname is not None:
            self.filename = fname

        try:
            self.qois = yaml.load(open(self.filename))
        except:
            raise

    def write(self,fname=None):
        """ write qoi configuration to yaml file

        Args:
            fname(str): file to write yaml from from.  If no argument is passed
               then use the filename attribute.  If the filename is set, then 
               the filename attribute is also set.
        """

        # set the attribute if not none
        if fname is not None:
            self.filename = fname

        # marshall attributes into a dictionary
        self.qoi_db = copy.deepcopy(self.qois)

        # dump to yaml file
        with open(self.filename,'w') as f:
            yaml.dump(self.qoi_db,f, default_flow_style=False)

class PotentialInformation(object):
    """ Read/Write Empirical Interatomic Potential Information

    pypospack uses yaml files to store configuration information for required
    for optimization routines.
    
    Attributes:
        filename(str): filename of the yaml file to read/write configuration
            file.  Default is pypospack.potential.yaml'
        elements(list): list of string of chemical symbols
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
