# -s*- coding: utf-8 -*-
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2017"
__license__ = "Simplified BSD License"
__version__ = "1.0"
"""This module contains the pyposmat engine for parameterization

The fitting engine is build in multiple class for managability of code 
and to divide the code into functionality.  In general, there are classes
which deal with reading/writing configuration files, classes to deal with
determination of simulations to be run from the quantities of interest.

Classes for Configuration Files
===============================
:obj:`pypospack.crystal.StructureDatabase`
:obj:`pypospack.qoi.QoiDatabase`
:obj:`pypospack.potential.PotentialInformation`

Classes for QOI Management
==========================
The definition of quantities of interest (QOI) are contained in a separate 
module where the definition of and calculation of specific quantities of 
interest are inherited from the :obj:`pypospack.qoi.QuantityOfInterest` 
prototype class.
:obj:`pypospack.qoi.QoiManager`
:obj:`pypospack.qoi.QoiDatabase`

Classes for LAMMPS Simulation Management
========================================
The LAMMPS simulation manager are defined as series of defined simulations,
which are run sequentially in the order required.  The simulations are defined
as tasks, which are combined togetheer to create a workflow.
:obj:`pypospack.lammps.SimulationManager`
:obj:`pypospack.tasks.lammps'

Classes for Structure Management
================================
:obj:`pypospack.crystal.StructureDatabase`
Classes for post-processing of data
===================================
"""

import time, yaml, copy
import os, shutil, subprocess
import importlib
import numpy as np
import scipy.stats

import pypospack.potential as potential
import pypospack.qoi as qoi
import pypospack.lammps as lammps
import pypospack.crystal as crystal

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


class PypospackFittingError(Exception):
    """Exceptional handling class for potential fitting"""
    def __init__(self,value):
        self.value = value

    def __self__(self):
        return repr(self.value)

class AbstractFittingEngine(object):
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
        structure_info(pypospack.potfit.StructureDatabase)
        potential_info(pypospack.potfit.PotentialInformation)
        qois_info(pypospack.potfit.QoiDatabase)
        qoi_manager(pypospack.qoi.QoiManager)
        simulation_manager(pypospack.lammps.SimulationManager)
    """

    def __init__(self,
            fname_config_structures = "pypospack.structure.yaml",
            fname_config_potential = "pypospack.potential.yaml",
            fname_config_qoi = "pypospack.qoi.yaml",
            fname_results = "results.out",
            fname_log = "pypospack.log",
            random_seed = None,restart = False):

        self.supported_qoi = list(get_supported_qois())
        self.supported_potentials = list(potential.get_supported_potentials())

        self.fname_config_structures = fname_config_structures
        self.fname_config_potential = fname_config_potential
        self.fname_config_qoi = fname_config_qoi
        self.fname_out_results = fname_results
        self.fname_out_log = fname_log

        self.obj_potential = None
        self.obj_lammps_sim_manager = None
        self.obj_qoi_manager = None

        self.restart = restart
        if self.restart is True:
            self.run_from_restart()

        # initialize output file
        self._initialize_log_file(fname_log)
        self._initialize_results_file(fname_results)

        # read in the configuration files
        self._process_input_file_structures(self.fname_config_structures)
        self._process_input_file_potential(self.fname_config_potential)
        self._process_input_file_qoi(self.fname_config_qoi)

        # check to see if structures required for qois exist
        self._check_qoi_structures_in_structure_db()

        # check external software
        self._get_lammps_binary()

        # determine random seed
        self._set_random_seed(random_seed)

        # configure obj_potential
        self._configure_potential()
        self._check_potential()

        #self._config_qoi_manager()
        self.qoi_manager = qoi.QoiManager(self.qoi_info)
        self.simulation_manager = lammps.SimulationManager()
        self.simulation_manager.add_required_simulations(
                required_simulations = self.qoi_manager.get_required_simulations())
        self.simulation_manager.structure_info = self.structure_info
        self.simulation_manager.potential_info = self.potential_info

    @property
    def qoi_names(self):
        """(:obj:`list` of :obj:`str`): a list of the qoi names"""
        return list(self.qoi_info.qois.keys())

    @property
    def parameter_names(self):
        """(:obj:`list` of :obj:`str`): a list of the parameter names"""
        return list(self.potential_info.parameter_names)

    @property
    def free_parameters(self):
        """(:obj:`list` of :obj:`str`): a list of free parameters"""
        return list(self.potential_info.free_parameters)

    def run_from_restart():
        raise NotImplementedError("Restart method not implemented")

    def evaluate_parameter_set(self,param_dict):
        self.simulation_manager.evaluate_parameters(param_dict)
        #self.qoi_manager.calculate_qois(results)

    def get_variables_from_simulation_manager(self):
        sim_mgr = self.simulation_manager
        assert isinstance(sim_mgr, lammps.SimulationManager)


    def _set_random_seed(self,seed = None):
        """Set the random seed in numpy.
        
        Since the seed might be randomly generated, we get the random_seed 
        value not from the parameter passed, but from numpy itself.

        Args:
            seed(int): the value of the seed.  The value must be an integer. 
               If the seed is set to None, the seed is generated automatically.

        Returns:
            int: returns the value of the seed.
        """
        if seed is None:
            self._log("auto_generate_seed:True")
       
        np.random.seed(seed)
        self.random_seed = np.random.get_state()[1][0]

        self._log("random_seed:{}".format(self.random_seed))

        return self.random_seed
    
    def _get_lammps_binary(self):
        try:
            self.lammps_bin = os.environ['LAMMPS_BIN']
        except KeyError as e:
            err_msg = "environment variable LAMMPS_BIN not set"
            self._log(err_msg)
            raise PypospackFittingError(err_msg)
        except:
            raise

        return self.lammps_bin

    def _initialize_results_file(self,fname_results = None):
        if fname_results is not None:
            self.fname_out_results = fname_results
        self._log("file_results -> {}".format(self.fname_out_results))
        self._f_results = open(self.fname_out_results,'w')

    def _initialize_log_file(self,fname_log = None):
        if fname_log is not None:
            self.fname_out_log = fname_log
        self._f_log = open(self.fname_out_log,'w')
        self._log("POTENTIAL OPTIMIZATION SOFTWARE")
        self._log("file_log -> {}".format(fname_log))

    def _log(self,msg):
        print(msg)
        self._f_log.write("{}\n".format(msg))

    def _process_input_file_structures(
            self,
            fname_config_structures = None):
        """ process structure file

        The structure database file is a yaml formatted file.  The creation of 
        this file can be done manually or using the class object
        :obj:`pypospack.potfit.StructureDatabase`

        Args:
            fname_config_structures(str): the filename of the structure 
                configuration file.  If no argument is passed, the class
                will use the attribute fname_config_structures.

        Raises:
            PypospackFittingError
        """

        if fname_config_structures is not None:
            self.fname_config_structures = fname_config_structures
        self._log("file_config_structure <- {}".format(self.fname_config_structures))
        self.structure_info = crystal.StructureDatabase()
        self.structure_info.read(self.fname_config_structures)
        structure_db_passed = self.structure_info.check()
        if structure_db_passed is not True:
            self._log(structure_db_passed)
            raise PypospackFittingError(structure_db_passed)

    def _process_input_file_potential(self,fname):
        """Process the potential definition
        
        The definition of the interatomic potential is defined a by a 
        formalism.  This class provides the information required to marshal
        and unmarshal information from a yaml file.  The creation of this
        file can be done manually or using the class object
        :obj:`pypospack.potfit.PotentialInformation`

        Args:
            fname_config_structures(str): the filename of the structure 
                configuration file.  If no argument is passed, the class
                will use the attribute fname_config_structures.

        """
        if fname is not None:
            self.fname_config_potential = fname
        self._log("file_config_potential <- {}".format(self.fname_config_potential))
        self.potential_info = potential.PotentialInformation()
        self.potential_info.read(self.fname_config_potential)
        self.potential_info.check() # sanity check

    def _process_input_file_qoi(self,fname):
        """Process the input file for quantities of interest (QOI)
        
        Quantities of interest are calculated from a variety of different
        simaultions.  This class provides the information required to marshall
        and unmaarshal information from a yaml file.  The creation of this
        file can be done manually or using the class object
        :obj:`pypospack.potfit.QoiDatabase`

        Args:
            fname_config_structures(str): the filename of the structure 
                configuration file.  If no argument is passed, the class
                will use the attribute fname_config_structures.

        """
        if fname is not None:
            self.fname_config_qoi = fname
        self._log("file_config_qoi <- {}".format(self.fname_config_qoi))
        self.qoi_info = qoi.QoiDatabase()
        self.qoi_info.read(self.fname_config_qoi)
        self.qoi_info.check() # sanity check

    def _check_qoi_structures_in_structure_db(self):
        """check qoi structures in the structure database

        The Quantities of Interest are dependent upon the calculation of 
        various structures.  These structure protypes need to exist in the
        fitting database.

        Returns:
            - **has_required_structures**(*bool*): True if the structure
              has all the structures in the fitting database.  False, otherwise.
            - **err_msg**(*str*): error message
        
        Raises:
            ValueError
        """

        has_required_structures = True # initialize
        
        required_structures = self.qoi_info.get_required_structures()
        missing_structures = [] # initialize
        for s in required_structures:
            if not self.structure_info.contains(structure = s):
                has_required_structures = False
                missing_structures.append(s)

        if not has_required_structures:
            # log and raise if there was a problem
            err_msg = "For the calcualtion of QOI's the following structures"
            err_msg += "are not contained in the structure database:\n"
            err_msg += "\n".join(missing_structures)
            self._log(err_msg)
            raise PypospackFittingError(err_msg)
        else:
            # there were no problem, returning true
            return True

    def _configure_potential(self):
        potential_type = self.potential_info.potential_type
        symbols = self.potential_info.symbols

        potential_map = potential.get_potential_map()
        module_name = potential_map[potential_type][0]
        class_name  = potential_map[potential_type][1]

        try:
            module = importlib.import_module(module_name)
            klass = getattr(module,class_name)
            self.obj_potential = klass(symbols)
        except:
            raise

    def _check_potential(self):
        if not set(self.obj_potential.parameter_names) == set(self.parameter_names):
            err_msg = "the potential parameter file does not have the parameter names as the formal definition\n"
            err_fmt = "{:^10}{:^10}\n"
            err_msg += err_fmt.format('config','definition')
            all_parameters = list(set(self.obj_potential.parameter_names + self.parameter_names))
            for p in all_parameters:
                cond1 = p in self.parameter_names
                cond2 = p in self.self.obj_potential.parameter_names
                err_msg += err_fmt.format(cond1,cond2)
            self._log(err_msg)
            raise PypospackFittingError(err_msg)

if False:
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

        def contains(self,structure):
            """ check to see that the structure in the structure database

            Args:
                structure(str):

            Returns:
                bool: True if structure is in structure database. False if the 
                    structure is not in the structure database
            """

            return structure in self.structures.keys()

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

        def check(self):
            """sanity checks for the fitting database
            
            This method checks the fitting database for the following errors: (1)
            the structure database exists, (2) files for the structure database
            exist

            Returns:
                str: returns a string if there is a problem
            Raises:
                ValueError: if there is a problem with the configuration
            """

            src_dir = self.directory
           
            # check to see if the source directory exists
            if not os.path.exists(src_dir):
                err_msg = "cannot find simulation directory\n"
                err_msg += "\tcurrent_working_directory:{}\n".format(os.getcwd())
                err_msg += "\tstructure_db_directory:{}\n".format(src_dir)
                return err_msg
            
            # check to see if the source directory is a directory
            if not os.path.isdir(src_dir):
                err_msg = "path exists, is not a directory\n"
                err_msg += "\tcurrent_working_directory:{}".format(os.getcwd())
                err_msg += "\tstructure_db_directory:{}\n".format(src_dir)
                return err_msg

            # check to see if files exist in the source directory
            files_exist = True
            msg = "structure files are missing:\n"
            for name, v in self.structures.items():
                filename = os.path.join(src_dir,v['filename'])
                if not os.path.isfile(filename):
                    files_exist = False
                    msg += "\t{}:{}\n".format(name,filename)

            if not files_exist:
                return msg
            else:
                return True

        def get_structure_dict(self,name):
            """

            Args:
                name(str): name of the structure
            """

            structure_db_dir = self.directory
            structure_filename = self.structures[name]['filename']
            structure_dict = {}
            structure_dict['name'] = name
            structure_dict['filename'] = os.path.join(
                    structure_db_dir,
                    structure_filename)

            return copy.deepcopy(structure_dict)

    class QoiDatabase(object):
        """ Qoi Database 
       
            Contains methods for managing quantities of interests for a variety 
            of tasks such as marshalling/unmarshalling objects to and from a
            yaml file.  Also includes sanity tests.

            Args:
                filename(str): If this variable is set, this class will attempt to
                    unmarhsall the contents of the file into the class.  If this 
                    variable is set to None, then the class wil be initialized.
                    Default is None.
            Attributes:
                filename(str): file to read/write to yaml file.  By default this
                    attribute is set to 'pypospack.qoi.yaml'.
                qois(dict): key is qoi name, value contains a variety of values
        """
            
        def __init__(self, filename = None):
            # initialize attributes
            self.filename = 'pypospack.qoi.yaml'
            self.qois = {}

            if filename is not None:
                self.read(filename)

        @property
        def qoi_names(self):
            return list(self.qois.keys())

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

        def check(self):
            """ perform sanity check on potential configuration
            
            does the following checks:
            1.    check to see if the potential type is supported.

            Raises:
                ValueError: if the potential type is unsupported
            # initialize variable, if a check fails the variable will be set to 
            # false
            """

            passed_all_checks = True

            for k,v in self.qois.items():
                qoi_type = v['qoi']
                if qoi_type not in get_supported_qois():
                    raise ValueError(\
                        "unsupported qoi: {}:{}".format(
                            k,qoi_type))

        def get_required_structures(self):
            """ get required structures """

            required_structures = []
            for qoi_name, qoi_info in self.qois.items():
                structures = qoi_info['structures']
                for s in structures:
                    if s not in required_structures:
                        required_structures.append(s)

            return required_structures

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
