# -*- coding: utf-8 -*-
"""This module contains the pyposmat engine for parameterization"""
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2016,2017"
__license__ = "Simplified BSD License"
__version__ = "1.0"

import time
import copy, shutil, os.path
import os, subprocess

import numpy as np
import scipy.stats

import pyflamestk.lammps as lammps
import pyflamestk.base as base
import pyflamestk.qoi as qoi
import pyflamestk.pareto as pareto

class PyPosmatError(Exception):
  """Exception handling class for pyposmat"""
  def __init__(self, value):
    self.value = value

  def __str__(self):
    return repr(self.value)

class PyPosmatEngine(object):

    def __init__(self,
                 fname_config_pyposmat = "pyposmat.config",
                 fname_config_potential = "pyposmat.potential",
                 fname_config_qoi = "pyposmat.qoi",
                 is_read = True,
                 is_restart = False,
                 random_seed = None):
        """default constructor

        this is the default constructor for the PyPosmatEngine is also 
        provides the base functionality and is designed so that it can be
        subclassed and extended for different types of implementation.

        Args:
            fname_config_pyposmat (str): the location of the pyposmat.config 
            fname_config_potential (str): the location of pyposmat.potential
            fname_config_qoi (str): the location of pyposmat.qoi
            is_read (bool): if set to True, then the configuration files will
                read
        """
        self._fname_config_pyposmat = fname_config_pyposmat
        self._fname_config_potential = fname_config_potential
        self._fname_config_qoi = fname_config_qoi

        self._is_read = is_read        # TODO: is this variable even used?
        self._is_restart = is_restart  # if set to true, 
                                       # try to recover previous simulations

        if (random_seed is not None) and (not is_restart):
            self._random_seed = random_seed
            np.random.seed(random_seed)

        self._config_pyposmat = None   # config info for pyposmat
        self._config_potential = None  # config info for potential
        self._config_qoi = None        # config info for qoi

        # read configuration files
        self._read_config_pyposmat(fname_config_pyposmat)
        self._read_config_potential(fname_config_potential)
        self._read_config_qoi(fname_config_qoi)

        # configure results file
        self._fname_out = self._config_pyposmat.fname_results_out
        self._f_results = open(self._fname_out,'w')

        # configure log file
        self._fname_log = self._config_pyposmat.fname_log
        self._f_log = open(self._fname_log, 'w')

        self._log('configuration files read')
        self._log('\tconfig_pyposmat -> {}'.format(self._fname_config_pyposmat))
        self._log('\tconfig_potential -> {}'.format(self._fname_config_potential))
        self._log('\tconfig_qoi -> {}'.format(self._fname_config_qoi))
        self._log('setting results file location:')
        self._log('\tresults_file -> {}'.format(self._fname_out))
        self._log('setting log file location:')
        self._log('\tlog_file -> {}'.format(self._fname_log))


        self._supported_qois = ['a0','a1','a2','a3',
                                'alpha','beta','gamma',
                                'c11','c12','c44',
                                'bulk_modulus',
                                'shear_modulus',
                                'defect_energy',
                                'surface_energy',
                                'total_energy']

        self._supported_potentials = ['buckingham',
                                      'eam']

        self._structure_db = None
        self._dir_structure_db = None
        self._validate_structure_database()

        self._dir_lammps_sim_db = None
        self._validate_lammps_simulation_templates()

        self._lmps_bin = None
        self._lmps_script = None
        self._validate_external_software()

        self._qoi_names = None
        self._qoi_info = None
        self._validate_quantities_of_interest()

        self._potential = None
        self._param_names = None
        self._param_dict = {}
        self._param_info = None
        self._configure_potential()

        self._var_names = None
        self._var_dict  = None
        self._qoi_manager = None
        self._qoi = None
        self._qoi_err = None
        self._configure_qoi_manager()

        self._lmps_sim_manager = None
        self._configure_lammps_simulation_manager()

        self._check_potential_parameters()
        self._sampler_type = None

        self._error_names = ["{}.err".format(q) for q in self._qoi_names]
        self._names = self._param_names + self._qoi_names + self._error_names
        self._types = len(self._param_names) * ['param'] \
                + len(self._qoi_names) * ['qoi'] \
                + len(self._error_names) * ['err']

    @property
    def structure_db(self):
        return self._structure_db

    @property
    def dir_structure_db(self):
        return self._dir_structure_db

    @property
    def sampler_type(self):
        return self._sampler_type

    @property
    def random_seed(self):
        return self._random_seed

    @sampler_type.setter
    def sampler_type(self, str_type):
        assert type(str_type),str
        self._sampler_type = str_type

    def _log(self,msg):
        print(msg)
        self._f_log.write(msg+'\n')

    def _restart_simulations(self):
        """
        this protected method implements the mechanism for restarting
        simulations
        """
        return start_sim_id

    def sample_parameter_space(self,
                               n_simulations,
                               fname_results='results.out',
                               sampler_type=None,
                               fname_results_in=None):
        """
        Parameters:
        n_simulations - number of simulations
        fname_results - filename of where to put simulation results
        sampler_type - supported_types: uniform, kde
        fname_results_in - required for kde
        """

        start_sim_id = 0          # initialized to normally start at 0
        write_header_files = True # initialized to normally writer headers
        self._results = None      # initialize results
        self._kde_kernel = None   # initialize

        if fname_results_in is not None:
            self._fname_results_in = fname_results_in

        if sampler_type is not None:
            self._sampler_type = sampler_type

        if self._is_restart is False:
            self._log("No restart, starting from sim_id:0")
            f = open(fname_results,'w')
        else:
            # restart requested, attempt to restart simulation
            self._log("Attempting simulation restart")
            try:
                n, t, self._results = self._read_results_file(fname_results)

                # get the last sim_id based on the last simulation result
                last_sim_id = self._results[len(self._results)-1][0]
                self._log("last sim_id:{}".format(last_sim_id))

                # set the simulaton id
                start_sim_id = last_sim_id + 1

                # append from the last file
                f = open(fname_results,'a')
                write_header_files = False
            except FileNotFoundError:
                # open the the file to write, since file was not found              
                f = open(fname_results,'w')
                self._log("simulation files not found, restarting from beginning")

        # header line strings
        header_line = ", ".join(['sim_id'] + self._names)
        types_line = ", ".join(['sim_id'] + self._types)

        # write headers
        if write_header_files:
            f.write(header_line + "\n")
            f.write(types_line + "\n")

        # do simulations
        for sim_id in range(start_sim_id,n_simulations):
            param_dict = None
            if self._sampler_type == 'uniform':
                param_dict = self.get_uniform_sample()
            elif self._sampler_type == 'kde':
                param_dict = self.get_kde_sample()
            else:
                raise PyPosmatError('unknown sampler type, \'{}\''.format(self._sampler_type))
            n,t,v = self.evaluate_parameter_set(param_dict)
            data = [sim_id] + v
            self._log("{}".format(sim_id))
            data_line = ", ".join(str(s) for s in data)

            f.write(data_line + '\n')
            if self._results is None:
                self._results = [data]
            else:
                self._results.append(data)
        self._results = np.array(self._results)
        f.close()

    def get_uniform_sample(self):

        # initialize param dict
        param_dict = {}
        for pn in self._param_names:
            param_dict[pn] = None

        for pn in self._param_names:
            if self._param_info[pn]['type'] == 'uniform':
                a = self._param_info[pn]['info'][0]
                b = self._param_info[pn]['info'][1]
                param_dict[pn] = np.random.uniform(a,b)
            elif self._param_info[pn]['type'] == 'static':
                param_dict[pn] = self._param_info[pn]['info'][0]

        # constrainted variables            
        for pn in self._param_names:
            if self._param_info[pn]['type'] == 'equals':
                info = self._param_info[pn]['info'][0]
                for p in self._param_names:
                    if p in info:
                        info = info.replace(p,"{}".format(param_dict[p]))
                param_dict[pn] = eval(info)

        return param_dict

    def get_free_parameter_list(self):
        self._free_param_list = []
        for pn in self._param_names:
            if self._param_info[pn]['type'] == 'uniform':
                self._free_param_list.append(pn)

        return self._free_param_list

    def initialize_kde_sampler(self,fname_in):
        f = open(fname_in,'r')
        lines = f.readlines()
        f.close()

        self._kde_names = lines[0].strip().split(',')
        self._kde_names = [str(v.strip()) for v in self._kde_names]

        self._kde_types = lines[1].strip().split(',')
        self._kde_types = [str(v).strip() for v in self._kde_types]

        datas = []
        for i in range(2,len(lines)):
            line = lines[i]
            line = line.strip().split(',')
            line = [float(v) for v in line]
            datas.append(line)
        datas = np.array(datas)

        # construct free parameter list
        free_param_list = self.get_free_parameter_list()
        self._kde_free_param_indx = []
        for i,v in enumerate(self._kde_names):
            if v in free_param_list:
                self._kde_free_param_indx.append(i)
        # DEBUG
        free_params = datas[:,self._kde_free_param_indx]
        self._kde_kernel = scipy.stats.gaussian_kde(free_params.transpose())

    def get_kde_sample(self):
        if self._kde_kernel is None:
            self.initialize_kde_sampler(self._fname_results_in)

        # initialize param dict
        param_dict = {}
        for pn in self._param_names:
            param_dict[pn] = None

        # construct free parameter list
        free_param_list = self.get_free_parameter_list()

        is_good = False
        while not is_good:
            # sample free parameters from kde
            free_params = self._kde_kernel.resample(size=1)
            for i,pn in enumerate(free_param_list):
                param_dict[pn] = free_params[i,0]

            # static variables
            for pn in self._param_names:
                if self._param_info[pn]['type'] == 'static':
                    param_dict[pn] = self._param_info[pn]['info'][0]

            # constrained variables
            for pn in self._param_names:
                if self._param_info[pn]['type'] == 'equals':
                    info = self._param_info[pn]['info'][0]
                    for p in self._param_names:
                        if p in info:
                            info = info.replace(p,"{}".format(param_dict[p]))
                    param_dict[pn] = eval(info)

            # check parameter constraints
            if param_dict['MgMg_rho'] < 0.:
                is_good = False
            elif param_dict['MgO_rho'] < 0.:
                is_good = False
            elif param_dict['OO_rho'] < 0.:
                is_good = False
            else:
                is_good =True

        return param_dict

    def evaluate_parameter_set(self, param_dict):
        #if not self._is_valid_parameter_set(param_dict):
        #    self._log("param_dict",param_dict)
        #    err_msg = "invalid parameter set
        #    raise PyPosmatError(err_msg) 

        # run the lammps simulations
        self._lmps_sim_manager.evaluate_parameter_set(param_dict)

        # extract pertinent variables from lammps simulations
        self._var_dict = self._lmps_sim_manager.variable_dict

        # calculate the qois
        self._qoi_manager.calculate_qois(self._var_dict)

        results_param = [param_dict[p] for p in self._param_names]
        results_qoi = [self._qoi_manager.qoi_values[q] for q in self._qoi_names]
        results_err = [self._qoi_manager.qoi_errors[q] for q in self._qoi_names]

        results = list(results_param) + list(results_qoi) + list(results_err)
        return self._names, self._types, results

    def _read_results_file(self, fname):
        """ reads the results file from a previous simulation
        
        Args:
            fname: the name of the filename to be read.

        Returns:
            A list of lists

        Raises:
            FileNotFoundError: the filename fname was not found.
        """
        assert type(fname),str

        write_header_files = False
        f = open(fname,'r')
        lines = f.readlines()
        f.close()

        # read the results file into memory
        results = None
        for i,line in enumerate(lines):
            if i in [0,1]:
                # process header lines
                lines[i] = [t.strip() for t in line.split(',')]
            else:
                # process simulation result lines
                lines[i] = [t.strip() for t in line.split(',')]
                lines[i] = [float(t) for t in lines[i]]
                lines[i][0] = int(lines[i][0]) # set sim_id

                if results is None:
                     # if there are no results, then add the first one
                    results = [[t for t in lines[i]]]

                else:
                    # add results
                    results.append([t for t in lines[i]])

        result_names = lines[0]
        result_types = lines[1]
        return result_names, result_types, results

    def _is_valid_parameter_set(self, param_dict):
        parameter_set_valid = True
        for pn in self._parameter_names:
            if pn not in self.param_dict.keys():
                parameter_set_valid = False
                self._log("ERROR: parameter set does not have {}".format(pn))
        return parameter_set_valid

    def create_lammps_simulations():
        raise NotImplementedError

    def _configure_potential(self):
        self._log('configuring the potential')
        if self._config_potential.potential_type == 'buckingham':
            self._log('\tbuckingham_potential')
            symbols = self._config_potential.symbols
            self._potential = lammps.BuckinghamPotential(symbols)
        elif self._config_potential.potential_type == 'tersoff':
            self._log('\ttersoff potential')
            symbols = self._config_potential.symbols
            self._potential = lammps.TersoffPotential(symbols)
        else:
            err_msg = 'potential_type is not supported. {}'
            err_msg = err_msg.format(self._config_potential.potential_type)
            raise ValueError(err_msg)
        self._potential.symbols = self._config_potential.symbols
        self._log('\tsymbols:' + ','.join(self._potential.symbols)) 
        self._param_names = self._potential.parameter_names
        self._param_info = self._config_potential.parameter_info
        self._log('\tparam_names:' + ','.join(self._param_names))
 
    def _read_config_pyposmat(self, fname):
        self._config_pyposmat  = PyPosmatConfigFile(fname_config=fname, 
                                                    is_read = True)
      
    def _read_config_potential(self, fname):
        self._config_potential = PotentialConfigFile(fname_config=fname, is_read = True)
       
    def _read_config_qoi(self, fname):
        self._config_qoi = QoiConfigFile(fname_config=fname, is_read = True)

    def _check_lammps_configuration(self):
        pass

    def _validate_structure_database(self):
        self._log('validating the structure database...')
        self._dir_structure_db = self._config_pyposmat.dir_structure_db
        if os.path.isdir(self._dir_structure_db):
            str_out = '\tdir_structure_db --> {}'
            str_out = str_out.format(self._dir_structure_db)
            self._log(str_out)
        else:
            err_msg = "could not find dir_structure_db[{}]"
            err_msg = err_msg.format(self._dir_structure_db)
            raise PyPosmatError(err_msg)

        # s_name = structure name
        # f_name = filename for the structure
        # s_type = filetype for the strucuture, vasp/lmmps
        for s_name in self._config_pyposmat.structure_names:
            f_name = self._config_pyposmat.structure_info[s_name][0]
            s_type = self._config_pyposmat.structure_info[s_name][1]

            full_path_name = os.path.join(self.dir_structure_db,f_name)

            # for lmmps type structures
            if s_type == 'lmmps':
                if os.path.isfile(full_path_name):
                    out_msg = '\tlmps: {} -> {}'
                    out_msg = out_msg.format(s_name,f_name)
                    self._log(out_msg)
                else:
                    err_msg = "file({}) for structure({}) does not exist"
                    err_msg = err_msg.format(f_name,s_namee)
                    raise PyPosmatError(err_msg)

            # for vasp type structures
            elif s_type == 'vasp':
                if os.path.isfile(full_path_name):
                    out_msg = '\tvasp: {} -> {}'
                    out_msg = out_msg.format(s_name,f_name)
                    self._log(out_msg)
                else:
                    err_msg = "file({}) for structure({}) does not exist"
                    err_msg = err_msg.format(f_name,s_name)
                    raise PyPosmatError(err_msg)
            else:
                # if structure type is not supported, raise an error
                # structure types suppored = 'lmmps'
                err_msg = "structure_type({}) is not supported"
                err_msg = err_msg.format(s_type)
                raise PyPosmatError(err_msg)

            self._structure_db = self._config_pyposmat.structure_info

    def _validate_lammps_simulation_templates(self):
        self._log('checking to see if lammps simulation database directory exists:')
        self._dir_lammps_sim_db = self._config_pyposmat.dir_lammps_sim_db
        if os.path.isdir(self._dir_lammps_sim_db):
            self._log('\tdir_lammps_sim_db -> {}'.format(self._dir_lammps_sim_db))
        else:
            err_msg = 'dir_lammps_sim_db[{}] does not exist'
            err_msg = err_msg.format(self._dir_lammps_sim_db)
            raise PyPosmatError(err_msg)

        # st_name = simulation template name
        # st_dir  = simulation template directory
        self._log('checking to see if all lammps simulation templates exist:')
        for st_name in self._config_pyposmat.lammps_sim_names:
            st_dir = self._config_pyposmat.lammps_sim_info[st_name]
            full_st_path = os.path.join(self._dir_lammps_sim_db, st_dir)
            if os.path.isdir(full_st_path):
                str_out = '\tlmmps_sim: {} -> {}'
                str_out = str_out.format(st_name,st_dir)
                self._log(str_out)
            else:
                err_msg = 'sim_template[{}] cannot find location[{}]'
                err_msg = err_msg.format(st_name,st_dir)
                raise PyPosmatError(err_msg)

    def _validate_external_software(self): pass

    def _validate_quantities_of_interest(self):
        self._log('validating quantities of interest.')
        self._qoi_info = self._config_qoi.qoi_info
        self._qoi_names = self._config_qoi.qoi_names
        # build a list of unique structures based upon qoi info
        structures = []
        for k in self._qoi_info.keys():
            new_structures = self._qoi_info[k]['structure']
            for s in new_structures:
                if s not in structures:
                    structures.append(s)

        # confirm that structures requested in the quantities of interest
        # are in the fitting database
        self._log('checking to see if qoi structures in structure_db.')
        for s in structures:
            if s in self._structure_db.keys():
                str_out = "\tqoi_structure: {} ... passed"
                str_out = str_out.format(s)
                self._log(str_out)
            else:
                err_out = "qoi_structure({}) not in in structure_db."
                err_out = err_out.format(s)
                raise PyPosmatError(err_out)
              
        # build list of qoi variables
        variables = []
        for qoi_key in self._qoi_info.keys():
            new_variable = self._qoi_info[qoi_key]['variable']
            if new_variable not in variables: 
                variables.append(new_variable)
              
        # confirm that the variables requested in the quantities of interest
        # can be calculated with the known simulation templates
        self._log('checking to confirm that vars required')
        for variable in variables:
            if variable in self._supported_qois:
                msg_out = "\tqoi_variable({}) is in qoi_db"
                msg_out = msg_out.format(variable)
                self._log(msg_out)
            else:
                err_out = "qoi_variable({}) is not in qoi_db"
                err_out = err_out.format(variable)
                raise PyPosmatError(err_out)

    def _configure_qoi_manager(self):
        self._log('setting up the qoi manager')
        self._qoi_manager = qoi.QoiManager()
        self._qoi_manager.qoi_definitions = self._qoi_info
        self._qoi_manager.process_qoi_definitions()
        self._var_names = self._qoi_manager.required_variable_names
        self._log('\tinitialized qoi manager')

    def _configure_lammps_simulation_manager(self):
        assert type(self._potential),lammps.Potential
        self._log('setting up the lammps simulation manager')
        self._lmps_sim_manager = lammps.SimulationManager()
        self._lmps_sim_manager.variable_names = self._var_names
        self._lmps_sim_manager.determine_simulations()
        self._lmps_sim_manager.potential = self._potential 
        self._lmps_sim_manager.structure_db = self._structure_db
        self._lmps_sim_manager.dir_structure_db = self._dir_structure_db
        self._lmps_sim_manager.dir_lmps_sim_db = self._dir_lammps_sim_db
        self._log('\tinitialized lammps simulation manager')

    def _check_potential_parameters(self): pass

class IterativeSampler(PyPosmatEngine):

    def __init__(self,
                 n_iterations,
                 n_simulations,
                 cull_type,
                 cull_param):

        # formatting strings
        self._iter_dir_format = "iter_{:03d}"
        self._fname_results_format = "results_{:03d}.out"
        self._fname_pareto_format = "pareto_{:03d}.out"
        self._fname_culled_format = "culled_{:03d}.out"

        # names of configuration files
        fname_config = 'pyposmat.config'
        fname_qoi = 'pyposmat.qoi'
        fname_potential = 'pyposmat.potential'

        PyPosmatEngine.__init__(self)

        # protected member variables
        self._n_iterations = n_iterations
        self._n_simulations = n_simulations
        self._cull_type = cull_type
        self._cull_param = cull_param

    @property
    def n_iterations(self):
        return self._n_iterations

    @n_iterations.setter
    def n_iterations(self, val):
        self._n_iterations = val

    @property
    def n_simulations(self):
        return self._n_simulations

    @n_simulations.setter
    def n_simulations(self, n_sims):
        self._n_simulations = n_sims

    def run(self,n_iterations=None):
        if n_iterations is not None:
            self._n_iterations = n_iterations

        for i_iter in range(self._n_iterations):
            self._log('starting iteration loop {}'.format(i_iter))
            fname_results_out = self._fname_results_format.format(i_iter)
            fname_pareto_out  = self._fname_pareto_format.format(i_iter)
            fname_culled_out  = self._fname_culled_format.format(i_iter)
            # generations
            if i_iter == 0:
                # use uniform sampling the first time around
                self.sample_parameter_space(n_simulations = self._n_simulations,
                                            fname_results = fname_results_out,
                                            sampler_type = 'uniform')
            else:
                # use the culled results from the previous iteration
                fname_results_in = self._fname_culled_format.format(i_iter-1)

                self.sample_parameter_space(n_simulations = self._n_simulations,
                                            fname_results = fname_results_out,
                                            sampler_type = 'kde',
                                            fname_results_in = fname_results_in)

            sim_results = pareto.SimulationResults()
            sim_results.read_simulation_results(fname_sims=fname_results_out)
            sim_results.calculate_pareto_set()
            sim_results.calculate_culled_set(self._cull_type,
                                             self._cull_param)
            sim_results.write_pareto_set(fname_pareto_out)
            sim_results.write_culled_set(fname_culled_out)

class FileParameterSampler(PyPosmatEngine):
    """This file samples parameters from a parameter file."""

    def get_qoi_header(self):
        """returns the header string for the output file

            Returns:
                str: space separated string of labels for output
        """
        str_out = ""
        str_out += "sim_id "
        for param in self.potential_parameter_list:
            str_out += "{} ".format(param)
        str_out += "| "
        for qoi_key in self.qoi_list:
            str_out += "{} ".format(qoi_key)
            str_out += "{}_abserr ".format(qoi_key)
            str_out += "{}_nabserr ".format(qoi_key)
            str_out += "{}_sqerr ".format(qoi_key)
            str_out += "{}_nsqerr ".format(qoi_key)
        str_out += "\n"
        return str_out

    def get_new_parameter_set(self,param_name_array,param_set):
        new_param_set = {}

        for param in self.potential_parameter_list:
            idx = param_name_array.index(param)
            param_value = param_set[idx] 
            new_param_set[param] = param_value
        return new_param_set

       
    def run(self,fname_params = 'params.in', fname_results = 'results.out'):

        # member variables initialization
        self.fname_params = fname_params
        self.fname_results = fname_results
        self.params_lines = []
        self.param_names  = []
        self.params = []

        # local variable initialziation
        n_simulations = 0
        n_sim_failures = 0
        start_time = time.time()
        end_time = time.time()

        #check to see if output file exists and delete
        if os.path.exists(fname_results):
            if os.path.isfile(fname_results):
                os.remove(fname_results)

        #check to see if file exists and read into memory
        if os.path.exists(self.fname_params):
            if os.path.isfile(self.fname_params):
                f = open(self.fname_params)
                self.params_lines = f.readlines()
                f.close()

        # get number of lines in file
        n_lines = len(self.params_lines)
        n_simulations = n_lines - 1

        # read header line
        self.param_names = self.params_lines[0].strip().split(' ')

        # read the parameter values from file into memory
        for idx in range(1,n_lines):
            self.params.append([float(num) for num in self.params_lines[idx].split()])

        # open output file      
        f = open(self.fname_results,'w')
        f.write(self.get_qoi_header())

        # do simulations loop
        for i_simulation in range(n_simulations):
            print('evaluating param_set_id:{}'.format(i_simulation))

            new_param_set = self.get_new_parameter_set(self.param_names,
                                                       self.params[i_simulation])

            self.evaluate_parameter_set(new_param_set)

            try:
                self.run_all_lammps_simulations()
                self.calculate_qoi()
                new_qoi_set = self.get_new_qoi_set()
         
                str_out = ""
                str_out = "{} ".format(i_simulation)
                for param in self.potential_parameter_list:
                    str_out += "{} ".format(new_param_set[param])
                str_out += "| "

                for qoi_result in new_qoi_set:
                  str_out += "{} ".format(qoi_result)
                str_out += "\n"
                f.write(str_out)

            except PyPosmatError:
                n_sim_failures += 1
                print("simulation failed, skipping parameter set") 
            # end simulations loop

        f.close()

        end_time = time.time()
        total_time = end_time - start_time
        print("\n\n")
        print("Simulations run: {}".format(n_simulations))
        print("Simulations failed: {}".format(n_sim_failures))
        print("Total time required for simulations: {} s".format(total_time))
        # end of run()

class SimulationResults(object):
    """this class processes the simulation results
  
    Args:
      n_simulations (int): number of simulations read from the output file
      qoi_type (str): supported qoi types are
          'abserr' - absolute error
    """

    def __init__(self):
        """default constructor"""
        self._supported_qoi_err_types = ['abserr', 'sqerr']

        # filenames    
        self.fname_log_file = "pyposmat.log"
        self.fname_sim_results = None
        self.fname_pareto = None
        self.fname_cull= None

        # initialize variables [ATTRIBUTES]
        self._qoi_err_type = 'abserr' # qoi error type
        self._qoi_err_req = None
        self._n_sims = None # numer of simulations
        self._names = None # names of the
        self._types = None
        self._param_names = [] # array of parameter names
        self._qoi_names = [] # array of qoi names
        self._err_names = []

        # results
        # set to None also indicates that calculations have not been done
        self._pareto_set_ids = None          # indexed with self._results
        self._results = None                 # numpy array of all simulation data
        self._pareto = None                  # numpy array of the pareto set
        self._cull = None                    # numpy array of the culled pareto set

        self.performance_requirements = {}

        # filename handles
        self._file_log = None
        self._log_format = None
        self._open_log_file()
 
    def __del__(self):
        self._close_log_file()

    @property
    def n_simulations(self): 
        return self._n_sims
  
    @n_simulations.setter
    def n_simulations(self, nsims):
        self._n_sims = nsims

    @property
    def qoi_ref(self):
        return self._qoi_ref

    @qoi_ref.setter
    def qoi_ref(self, dict_qoi):
        assert type(dict_qoi),dict
        self._qoi_ref = dict_qoi

    @property
    def qoi_err_req(self):
        return self._perf_req

    @property
    def names(self):
        """List of str: contains a list of all string in the column"""
        return self._names

    @property
    def types(self):
        """List of str: contains either 'param','qoi','err'"""
        return self._types

    @property
    def qoi_err_type(self):
        return self._qoi_err_type

    @qoi_err_type.setter
    def qoi_err_type(self, qe_type):
        assert type(qe_type),str

        # check to see if qoi_err_type is supported
        if qe_type not in self._supported_qoi_err_types:
            err_msg = 'unsupported qoi error type: {}'
            err_msg = err_msg.format(qe_type)
            raise ValueError(err_msg)

        self._qoi_err_type = qe_type

    @property
    def qoi_names(self):
        return self._qoi_names

    @property
    def err_names(self):
        return self._err_names

    @qoi_names.setter
    def qoi_names(self, qnames):
        self._qoi_names = qnames

    @property
    def parameter_names(self):
        return self._param_names

    @parameter_names.setter
    def parameter_names(self, pnames):
        self._param_names = pnames
    
    @property
    def results(self):
        """numpy.array: numpy array of results"""
        return self._results

    @property
    def pareto(self): 
        """numpy.array: numpy array of pareto results"""
        return self._pareto
 
    @property
    def culled(self): 
        """numpy.array: numpy array of pareto results"""
        return self._culled

    # <--- some functions here to read configuration files
    def read_configuration_files(self,
                                 pyposmat_config = 'pyposmat.config',
                                 pyposmat_potential = 'pyposmat.potential',
                                 pyposmat_qoi = 'pyposmat.qoi'):
        """" reads the pyposmat configuration files
        
        Args:
            pyposmat_config (str,optional): defaults to pyposmat.config
            pyposmat_potential (str,optional): defaults to pyposmat.potential
            pyposmat_qoi (str,optional): defaults to pyposmat.qoi
      
        """
        self._read_config_pyposmat(pyposmat_config)
        self._read_config_potential(pyposmat_potential)
        self._read_config_qoi(pyposmat_qoi)

    def _read_config_pyposmat(self,fname):
        self._config_pyposmat = PyPosmatConfigFile(fname, True)

    def _read_config_potential(self,fname):
        self._config_potential = PotentialConfigFile(fname, True)

    def _read_config_qoi(self,fname):
        self._config_qoi = QoiConfigFile(fname, True)
        self._qoi_ref = self._config_qoi.qoi_ref_values

    # SOME FUNCTIONS HERE TO DEAL WITH APPLICATION LOGGING.          
    def _open_log_file(self):
        if self._file_log is None:
            self._file_log = open(self.fname_log_file,'w')
  
    def _close_log_file(self):
        self._file_log.close()
  
    def _log(self,msg):
        if self._log_format is None:
            self._file_log.write(msg + "\n")
            print(msg)
        else:
            msg = self._log_format.format(msg + "\n")
            print(msg)
    
    def write_pareto_set(self,fname_out='pareto.out'):
        """Write the pareto set to file.

        This function prints the calculated pareto set to file.

        Args:
            fname_out(str) - the filename (default: pareto.out)
        """

        # create header
        str_names = ", ".join(self._names) + "\n"
        str_types = ", ".join(self._types) + "\n"

        # create body
        str_body = ""
        for sim_result in self._pareto:
            str_body += ", ".join([str(num) for num in sim_result]) + "\n"
          
        # write results
        f = open(fname_out,'w')
        f.write(str_names)
        f.write(str_types)
        f.write(str_body)
        f.close()

    def write_culled_set(self,fname_out='culled.out'):
        # create header
        str_names = ", ".join(self._names) + "\n"
        str_types = ", ".join(self._types) + "\n"

        # create body
        str_body = ""
        for sim_result in self._culled:
            str_body += ", ".join([str(num) for num in sim_result]) + "\n"
          
        # write results
        f = open(fname_out,'w')
        f.write(str_names)
        f.write(str_types)
        f.write(str_body)
        f.close()
    
    def write_analysis_files(self,
                             dir_name = None,
                             fname_pareto = 'pareto.dat',
                             fname_culled = 'culled.dat',
                             is_write_pareto = True, 
                             is_write_culled_set = True):
        """
        writes a variety of analysis files
      
        Arguments:
      
        dir_name (str) - destination directory name in which to put files
        """

        if not (dir_name == None):
            self.working_path = dir_name
        else:
            # self.working_path stays the same
            pass

        # create directory if directory does not exist
        os.makedirs(dir_name, exist_ok=True)
        msg = "working path: {}".format(self.working_path)
        self.__log(msg)

        # write results of the pareto set
        if is_write_pareto == True:
            fname = os.path.join(dir_name,fname_pareto)
            self.__log("writing pareto set to {}".format(fname))
            self.__write_pareto_set(fname)
          
        # write results of the culled pareto set
        if is_write_culled_set == True:
            fname = os.path.join(dir_name,fname_pareto)
            self.__log("writing culled pareto set to {}".format(fname))
            self.__write_culled_set()

    def __read_file(self, fname, file_type):
 
        # read file into memory
        f_in = open(fname,'r')
        lines_in = f_in.readlines()
        f_in.close()

        if file_type == 'results':
            # read header lines
            self._names = [n.strip() for n in lines_in[0].strip().split(',')]
            self._types = [t.strip() for t in lines_in[1].strip().split(',')]
        elif file_type == 'pareto':
            # check to see if pareto header line is the same as the pareto line
            if self._names == [n.strip() for n in lines_in[0].strip().split(',')]:
                pass
            else:
                if self._names is None:
                    errmsg = "The results file must be read before the pareto file"
                    raise RuntimeError(errmsg)
                else:
                    errmsg = "The pareto names header does not match results file"
                    raise RuntimeError(errmsg)

            # check to see if pareto types header line is the same as the pareto line
            if self._types == [t.strip() for t in lines_in[1].strip().split(',')]:
                pass
            else:
                if self._types is None:
                    errmsg = "The results file must be read before the pareto file"
                    raise RuntimeError(errmsg)
                else:
                    errmsg = "the pareto types header does not match results file"
                    raise RuntimeError(errmsg)

        results = []
        for i in range(2,len(lines_in)):
            result =  [v.strip() for v in lines_in[i].strip().split(',')]
            for j,t in enumerate(self._types):
                if t == "sim_id":
                    result[j] = int(float(result[j]))
                else:
                    # everything else is a float
                    result[j] = float(result[j])
            results.append(result)

        # convert into numpy file
        if file_type == 'results':
            self._param_names = [self._names[i] for i,v in enumerate(self._types) if v == 'param']
            self._qoi_names = [self._names[i] for i,v in enumerate(self._types) if v == 'qoi']
            self._err_names = [self._names[i] for i,v in enumerate(self._types) if v == 'err']
            self._results = np.array(results)
        elif file_type == 'pareto':
            self._pareto = np.array(results)
        elif file_type == 'culled':
            self._culled = np.array(results)

    def read_simulation_results(self,
                                fname_sims, 
                                fname_pareto = None, 
                                fname_cull = None):
        """
        read simulations results from a file into a memory.
      
        Args:
            fname_sims (str): the filename containing the simulation results from
                LAMMPS simulations
            fname_pareto (str): the filename containing the pareto set results
            fname_cull (str): the filename contain culled pareto set results
        """
      
        self.fname_sims = fname_sims
        self.__read_file(fname_sims, 'results')
        
        # remove rows that have NaN as result
        rows_to_remove = []
        for i in range(1,self._results.shape[0]):
            if np.isnan(self._results[i,:]).any():
                rows_to_remove.append(i)
        self._results = np.delete(self._results,rows_to_remove,axis=0)
      
        if fname_pareto is not None:
            self.fname_pareto = fname_pareto
            self.__read_file(fname_pareto, 'pareto')
          
        if fname_cull is not None:
            self.fname_cull = fname_cull
            self.__read_file(fname_cull, 'culled')
          
    def _create_dataset_for_pareto_analysis(self, err_names=None):
        """iCreates a dataset for pareto analysis

        This method creates a dataset necessary for pareto analysis

        Arguments:
        err_names (list of str): - contains the identifiers for the error
        """

        print("creating dataset for pareto analysis")
         
        if err_names is None:
            err_names = self._err_names
            
        # get indices of error names
        err_idx = [self._names.index(n) for n in err_names]

        # select the sim_id column and err_names columns
        results_err = self._results[:,[0] + err_idx]
        results_abs_err = np.abs(results_err)

        # make dataset 
        n_row, n_col = results_abs_err.shape
        self._pareto_dataset = []
        for i_row in range(n_row):
            self._pareto_dataset.append(Datapoint(i_row))
            for i_col in range(n_col):
                number = results_abs_err[i_row,i_col]
                self._pareto_dataset[i_row].addNumber(-number)
 
    def calculate_pareto_set(self):

        self._create_dataset_for_pareto_analysis(err_names=self._err_names)
        bruteforce_algo(self._pareto_dataset)

        # mark pareto set
        pareto_set = []
        pareto_set_ids = []
        for s in self._pareto_dataset:
            if s.paretoStatus == 1:
                pareto_set_ids.append(s.id)
                pareto_set.append(s.vec)
          
        #pareto_set = -np.array(pareto_set)
        self._pareto = self._results[pareto_set_ids,:]

    def calculate_parameter_estimates(self,param_list):
        params = copy.deepcopy(param_list)
        self.param_estimates = {}
        for param in params:
            self.param_estimates[param] = {}
            self.param_estimates[param]['all'] = {}
            self.param_estimates[param]['pareto'] = {}
            self.param_estimates[param]['pareto_cull'] = {}

    def calculate_qoi_estimates(self,qoi_keys):
        qois = copy.deepcopy(qoi_keys)
        set_types = ['all','pareto','pareto_cull']
        self.qoi_estimates = {}
        for qoi in qois:
            self.qoi_estimates[qoi] = {}
            for set_type in set_types:
                #TODO
                mean = 0
                std  = 1
                self.qoi_estimates[qoi][set_type] = {}
                self.qoi_estimates[qoi][set_type]['mean'] = mean
                self.qoi_estimates[qoi][set_type]['std'] = std
          
    #--------------------------------------------------------------------------
    # methods for calculating the culled pareto set
    #--------------------------------------------------------------------------
    def calculate_culled_set(self,
                             cull_type="percentile",pct=80.,
                             qoi_err_threshold = None):
        """calculates the culled set
        
        This method applies applies performance thresholds upon the simulation.
        If qoi_error_threshold is set, the parameters cull_type and pct are 
        ignored.

        Args:
            cull_type (str): supports the different culling of the pareto set by 
                different mechanisms.  The current mechanisms are 'percentile'
                and 'pct_error'
            pct (float): is a float variable.
            qoi_err_threshold: the error threshold acceptable.

        Returns:
            Nothing

        Raises:
            RuntimeError: If any key in qoierr_threshold is not contained
                in the attribute error_names, it will check to see if
                the key value is contained in qoi_names and attempt to 
                change the key value.  If this is not successful, then
                a RuntimeError will be returned.
        """

        if qoi_err_threshold is not None:
             for k,v in qoi_err_threshold.items():
                 if k.endsin('.err'):
                     if k in self.error_names:
                         self.add_performance_constraint(k,v)
                     else:
                         raise RuntimeError('unknown performance constraint')
                 else:
                     if k in self.qoi_names:
                         new_k = "{}.err".format(k)
                         self.add_performance_constraint(new_k,v)
                     else:
                         raise RuntimeError('unknown performance constraint')
             self.apply_performance_constraints()
             return

        if cull_type == "percentile":
            self._calculate_culled_set_by_percentile(pct)
        elif cull_type == "pct_error":
            self._calculate_culled_set_by_percent_error(pct)
        else:
            raise RuntimeError("unknown cull_type")

    def _calculate_culled_set_by_percentile(self,pct_kept=80.):
        """
        Arguments:
        pct_kept (float, 10.0) - number between 1 and 100 indicating the 
        pct of simulations within the Pareto set which should be kept
                        
        Returns:
        
        a numpy array with observations indexed in rows, and parameters
        and quantities of interst indexed in columns.  The column index is the
        same as the array in "self.all_names".
        """

        # TODO:
        # A Newton-Ralphson method to get more accurate performance requirements
        #to prevent over culling of the Pareto set.
       
        if not(0 <= pct_kept <= 100.):
            errmsg = "pct_kept must be between 1 and 100, the value {} was passed."
            errmsg = errmsg.format(pct_kept)
            raise ValueError(errmsg)
        else:
            self.pct_kept = pct_kept
              
        err_keys = self._err_names
        self._perf_req = {}
        for err_key in err_keys:
            self._perf_req[err_key] = 0.
        n_sims, n_qoi = self._pareto.shape        
    
        # intialize variables
        pctl_threshold = 100        # searching for 100% within the Pareto set
                                    # to 0% in the pareto set
        is_culled = False           # intialize
        
        while not is_culled:
            rows_to_delete = []
            pctl_threshold -= 0.1
            # calculate percentile cutoffs
            for err_key in self._perf_req.keys():
                if pctl_threshold < 0:
                    errmsg = "While searching for the pctl_threshold, the \
                              percentile error dropped below zero resulting \
                              in an error."
                    raise ValueError(errmsg)
                else:
                    qoi_data = self.get_data(err_key, 'pareto','abserror')
                    cutoff = np.percentile(qoi_data,pctl_threshold)
                    self._perf_req[err_key] = cutoff

            # cull the pareto set by the performance requirements
            for err_key in self._perf_req.keys():        
                pareto = np.copy(self.pareto)
                for idx in range(n_sims):
                    ps = pareto[idx,:]

                    # determine if row needs to be deleted
                    is_delete_row = False
                    for qoi_name in self._perf_req.keys():
                        qoi_idx = self._names.index(qoi_name)
                        if ps[qoi_idx] > self._perf_req[qoi_name]:
                            is_delete_row = True 
                            
                    # add row for deletion if requirements met.
                    if is_delete_row:
                        rows_to_delete.append(idx)
            
            # check to see if the pareto set has been sufficiently culled
            n_culled = len(rows_to_delete)
            pct_culled = float(n_culled)/float(n_sims)
            if pct_kept/100. > 1 - pct_culled:
                is_culled = True

        self._culled = np.delete(self._pareto,
                                 rows_to_delete,
                                 axis=0)                
        return self._culled.copy()
            
    def _calculate_culled_set_by_percent_error(self,pct_kept=80.):
        """
        Keeps members of the pareto set, which are (1+pct) over the
        reference value
        
        Arguments:
            pct (float, 8.0) - number indicating the cutoff within which
                members of the Pareto set should be kept.
        """
        if (pct_kept) <= 0:
            errmsg = "pct_kept must be between 1 and 100, the value {} was passed."
            errmsg = errmsg.format(pct_kept)
            raise ValueError(errmsg)
        
        # calculate performance constraints.
        pct_kept = float(pct_kept)             # force casting into float
        pct_kept = pct_kept/100.
        self._perf_req = {k:pct_kept*v for k,v in self._qoi_ref.items()}
        
        self.apply_performance_constraints()

    def _calculate_culled_set_by_performance_requirements(self, perf_req=None):
        if perf_req is None:
            perf_req = self._perf_req
        else:
            assert type(perf_req), dict

        self.apply_performance_constraints()

    def add_performance_constraint(self,metric_name,metric_value):
        assert type(metric_name),str
        assert type(metric_value),float

        if self._perf_req is None:
            self._perf_req = {}
        self._perf_req[metric_name] = metric_value

    def apply_performance_constraints(self):
        # start with the full pareto set and then remove elements which do 
        # not meat the performane criteria
        n_sims, n_qoi = self._pareto.shape
        self._culled = np.copy(self._pareto)

        #determine which rows to delete
        rows_to_delete = []
        for idx in range(n_sims):
            ps = self._culled[idx,:]
            is_delete_row = False
            for qoi_name in self._perf_req.keys():
                err_name = qoi_name + ".err"
                err_idx = self._names.index(err_name)
                if np.abs(ps[err_idx]) > self._perf_req[qoi_name]:
                    is_delete_row = True 
            if is_delete_row:
                rows_to_delete.append(idx)
        
        # remove rows which do not meet performance criteria
        self._culled = np.delete(self._culled,
                                 rows_to_delete,
                                 axis=0)

    def get_data(self, name, ds_type, err_type = 'abserr'):
        """get data

        Args:
            name (str): string of parameter or quantity of interest
            ds_type (str): string of which dataset we are taking the data from.
                The ds_types which are supported are: all, pareto, pareto_culled
          
        Returns:
            (numpy.array): a numpy array of the data asked for
        """
        idx  = self._names.index(name)

        # get data by dataset
        data = None # initialize
        if ds_type == 'results':
            data = self._results[:,idx]
        elif ds_type == 'pareto':
            data = self._pareto[:,idx]
        elif ds_type == 'culled':
            data = self._culled[:,idx]

        if self._types[idx] == 'err':
            # transform errors if necessary
            if err_type == 'err':
                # no transformation required
                data = self._results[:,idx]
            elif err_type == 'abserr':
                # transform for absolute errors
                data = np.abs(self._results[:,idx])
        else:
            # tranformation not necessary
            data = self._results[:,idx]

        return copy.deepcopy(data)
   
    def create_all_pareto_plots(self,qoi_list):
        for i, qoi_name_i in enumerate(qoi_list):         
            for j, qoi_name_j in enumerate(qoi_list):
                if i < j and i != j:
                    print("{} {}".format(qoi_name_i,qoi_name_j))
                    x_label = qoi_name_i
                    y_label = qoi_name_j
                    x_idx = self.all_names.index(x_label)
                    y_idx = self.all_names.index(y_label)
                    pareto_front = pareto_frontier_2d(self.np_all_sims[:,x_idx],
                                          self.np_all_sims[:,y_idx],
                                          maxX = False, maxY= False)
                    fig, ax = plt.subplots()
                    ax.scatter(self.np_all_sims[:,x_idx],
                           self.np_all_sims[:,y_idx],
                           label='dominated')
                    ax.scatter(self.pareto_set[:,self.qois.index(x_label)],
                           self.pareto_set[:,self.qois.index(y_label)],
                           label = 'pareto', color='y',)
                    ax.plot(pareto_front[0],
                            pareto_front[1],
                            color='r',
                            linewidth=2)
                    legend = ax.legend(loc="upper right")
                    plt.axis([min(pareto_front[0]),
                              max(pareto_front[0]),
                              min(pareto_front[1]),
                              max(pareto_front[1])])
                    plt.xlabel(x_label)
                    plt.ylabel(y_label)
                    plt.show()

class PyPosmatConfigFile:
    def __init__(self,
                 fname_config = "pyposmat.config",
                 is_read = True):
        self._pyposmat_out_fname = "results.out"
        self._pyposmat_log_fname = "pyposmat.log"
        self._fname = fname_config

        self._structure_names = []
        self._structure_info = None
        self._dir_structure_db = None

        self._lmps_sim_names = []
        self._lmps_sim_info = None
        self._dir_lammps_sim_db = None

        self._lmps_bin = None
        self._lmps_exe_script = None

        if is_read is True:
            self.read()

    @property
    def fname_results_out(self):
        return self._pyposmat_out_fname

    @property
    def fname_log(self):
        return self._pyposmat_log_fname

    @property
    def structure_names(self):
        return self._structure_names

    @property
    def structure_info(self):
        return self._structure_info

    @property
    def lammps_bin(self):
        return self._lmps_bin

    @property
    def lammps_exe_script(self):
        return self._lmps_exe_script

    @property
    def dir_structure_db(self):
        return self._dir_structure_db

    @property
    def lammps_sim_names(self):
        return self._lmps_sim_names

    @property
    def dir_lammps_sim_db(self):
        return self._dir_lammps_sim_db

    @property
    def lammps_sim_info(self):
        return self._lmps_sim_info

    def read(self):
        f = open(self._fname,'r')
        lines = f.readlines()
        f.close()
        for line in lines:
            line = line.strip()
            # skip if line is a comment
            # skip if line is empty or whitespace
            if line.startswith('#') or line == '':
                pass
            else:
                config_info = line.split('=')
                config_var = config_info[0].strip()
                config_param = config_info[1]
                if config_var == 'structure':
                    self._add_structure(config_param)
                elif config_var == 'lmps_sim_type':
                    self._add_lammps_simulation_type(config_param)
                elif config_var == 'pyposmat_out_fname':
                    self._pyposmat_out_fname = config_param.strip()
                elif config_var == 'pyposmat_log_fname':
                    self._pyposmat_log_fname = config_param.strip()
                elif config_var == 'pyposmat_type':
                    self._pyposmat_type = config_param.strip()
                elif config_var == 'lmps_bin':
                    self._lmps_bin = config_param.strip()
                elif config_var == 'lmps_exe_script':
                    self._lmps_exe_script = config_param.strip()
                elif config_var == 'dir_structure_db':
                    self._dir_structure_db = config_param.strip()
                elif config_var == 'dir_lammps_sim_db':
                    self._dir_lammps_sim_db = config_param.strip()
                else:
                    err_msg = "Unknown configuration command, {}"
                    err_msg = err_msg.format(config_var)
                    raise PyPosmatError(err_msg)

    def _add_structure(self,config_param):
        s_name  = config_param.split(',')[0].strip()
        s_fname = config_param.split(',')[1].strip()
        s_type  = config_param.split(',')[2].strip()
        self._structure_names.append(s_name)
        if self._structure_info is None:
            self._structure_info = {}
        self._structure_info[s_name] = [s_fname,s_type]

    def _add_lammps_simulation_type(self,config_param):
        sim_name = config_param.split(',')[0].strip()
        sim_location = config_param.split(',')[1].strip()
        self._lmps_sim_names.append(sim_name)
        if self._lmps_sim_info is None:
            self._lmps_sim_info = {}
        self._lmps_sim_info[sim_name] = sim_location


class ConfigFile:
  def __init__(self, fname_config = "pyposmat.config",is_read = True):
    self.fname_config = fname_config
    if is_read == True:
      self.read()
      
  def get_pyposmat_info(self):
      self.pyposmat_type = self.config_dict['pyposmat_type']
      self.pyposmat_log  = self.config_dict['pyposmat_log_fname']
      self.pyposmat_out  = self.config_dict['pyposmat_out_fname']
      return {'pyposmat_type': self.pyposmat_type,
              'pyposmat_log' : self.pyposmat_log,
              'pyposmat_out' : self.pyposmat_out}

  def get_number_of_objective_functions(self):
    n_qoi         = len(self.config_dict['qoi'])
    n_qoi_weights = len(self.config_dict['qoi_weights'])
    n_qoi_targets = len(self.config_dict['qoi_target'])
    
    if (n_qoi == n_qoi_targets) and (n_qoi == n_qoi_targets):
        # check to see if qoi, weights, and values are equal
        return n_qoi
    else:
        err_msg = 'n_qoi({}), n_qoi_weights({}), n_qoi_targets({}) are different.'
        err_msg.format(n_qoi, n_qoi_weights, n_qoi_targets)
        raise ValueError(err_msg)
      
  def getNumberOfObjectiveFunctions(self):
    n_qoi         = len(self.config_dict['qoi'])
    n_qoi_weights = len(self.config_dict['qoi_weights'])
    n_qoi_targets = len(self.config_dict['qoi_target'])
    return n_qoi

  def getDakotaKwargs(self):
    dakota_kwargs = {'dakota_bin'      :self.config_dict["dakota_bin"],
                     'dakota_fname_in' :self.config_dict["dakota_fname_in"],
                     'dakota_fname_out':self.config_dict["dakota_fname_out"],
                     'dakota_fname_err':self.config_dict["dakota_fname_err"]}
    return dakota_kwargs

  def read(self):
    self.config_dict = read_pyposmat_config_file(self.fname_config)
    return copy.deepcopy(self.config_dict)
    
  
def read_pyposmat_config_file(fname_config):
  config_dict = {}
  config_dict['structures'] = {}
  config_dict['lmps_sim_type'] = {}
  config_dict['qoi_list'] = []
  config_dict['qoi'] = {}
  config_dict['qoi_weights'] = {}
  config_dict['qoi_targets'] = {}
  config_dict['qoi_normfactor'] = {}
  file = open(fname_config,'r')
  for idx, line in enumerate(file.readlines()):
    if line.strip().startswith('#'):
      # skip if the line is a comment
      pass
    elif line.strip() == '':
      # skip if the line is empty or only contains white space
      pass
    else:
      config_info = line.split('=')
      if len(config_info) > 0:
        config_var = config_info[0].strip()
        config_param = config_info[1].strip().strip('"').strip()
        if config_var == 'structure':
          # structure = structure_name, structure_fname, structure_type
          structure_name  = config_param.split(',')[0].strip()
          structure_fname = config_param.split(',')[1].strip()
          structure_type  = config_param.split(',')[2].strip()
          config_dict['structures'][structure_name] = {}
          config_dict['structures'][structure_name]['name'] = structure_fname
          config_dict['structures'][structure_name]['type'] = structure_type
        elif config_var == "lmps_sim_type":
          sim_type     = config_param.split(',')[0].strip()
          sim_location = config_param.split(',')[1].strip()
          config_dict['lmps_sim_type'][sim_type] = sim_location
        elif config_var == "qoi_list":
          config_dict["qoi_list"] = [word.strip() for word in config_param.split(',')]
        elif config_var == "qoi":
          args = config_param.split(',')
          n_structures = len(args) - 2
          qoi_name = args[0].strip()
          qoi_type = args[1].strip()
          config_dict['qoi'][qoi_name] = {}
          config_dict['qoi'][qoi_name]['type'] = qoi_type
          config_dict['qoi'][qoi_name]['structures'] = [args[i].strip() for i in range(2,2 + n_structures)]
        elif config_var == "qoi_weights":
          qoi_name   = config_param.split(',')[0].strip()
          qoi_weight = float(config_param.split(',')[1].strip())
          config_dict['qoi_weights'][qoi_name] = qoi_weight
        elif config_var == "qoi_target":
          qoi_name       = config_param.split(',')[0].strip()
          qoi_target_val = float(config_param.split(',')[1].strip())
          config_dict['qoi_targets'][qoi_name] = qoi_target_val
        elif config_var == "qoi_normfactor":
          qoi_name       = config_param.split(',')[0].strip()
          qoi_normf      = float(config_param.split(',')[1].strip())
          config_dict['qoi_normfactor'][qoi_name] = qoi_normf
        else:
          config_dict[config_var] = config_param
  return config_dict

class PotentialConfigFile:
    def __init__(self,
                 fname_config = "pyposmat.potential",
                 is_read = True):

        # class member variables
        self.fname_config = fname_config
        self.config_dict  = {}
        self._elements = {}
        self._symbol_list = None
        self._param_names = []
        self._param_info = {}
        self._pair_type = {}
        self.pair = {}
    
        if is_read == True:
            self.read()

    @property
    def symbols(self):
        return self._symbol_list

    @property
    def elements(self):
        return self._elements

    @property
    def potential_type(self):
        # NO SUPPORT FOR HYBRID POTENTIALS
        # ALL PAIRS MUST BE OF THE SAME TYPE
        pair_types = []
        for k in self._pair_type.keys():
            pair_type = self._pair_type[k]['type']
            assert type(pair_type),str
            if pair_type not in pair_types:
                pair_types.append(pair_type)

        # there should be only one pair_type
        assert len(pair_types),1
        return pair_types[0]

    @property
    def parameter_info(self):
        return self._param_info

    def read(self):

        file = open(self.fname_config)
        lines = file.readlines()
        file.close()

        for line in lines:
            line = line.strip()
            if line.startswith('#') or line == '':
                pass
            else:
                keyword = line.split('=')[0].strip()
                params  = line.split('=')[1].strip().split()
                if keyword == 'potential_elements':
                    self._symbol_list = [p.strip() for p in params]
                    for param in params:
                        self._elements[param.strip()] = {}
                elif keyword == 'potential_charge':
                    self._set_potential_charge_param(params)
                elif keyword == 'potential_pair_type':
                    self._set_potential_pair_type(params)
                elif keyword == 'potential_3body_type':
                    self._set_potential_3body_type(params)
                elif keyword == 'potential_pair_param':
                    self._set_potential_pair_param(params)
                elif keyword == 'potential_3body_type':
                    self._set_potential_3body_param(params)

    def _set_potential_charge_param(self,params):
        params = [p.strip() for p in params]
        p_name = 'chrg_{}'.format(params[0])
        p_type = params[1]
        p_info = [params[i] for i in range(2,len(params))]

        self._set_param(p_name,p_type,p_info)

    def _set_potential_pair_param(self,params):
        params = [p.strip() for p in params]
        p_name = '{}{}_{}'.format(params[0],
                                  params[1],
                                  params[2])
        p_type = params[3]
        p_info = [params[i] for i in range(4,len(params))]

        self._set_param(p_name,p_type,p_info)

    def _set_potential_3body_param(self,params):
        params = [p.strip() for p in params]
        element1 = params[0]
        element2 = params[1]
        element3 = params[2]
        param_name = params[3]
        p_name = '{}{}{}_{}'.format(element1,
                                    element2,
                                    element3,
                                    param_name)
        p_type = params[4]
        p_info = [params[i] for i in range(5,len(params))]
        self._set_param(p_name,p_type,p_info)

    def _set_param(self,p_name,p_type,p_info):
        self._param_names.append(p_name)
        self._param_info[p_name] = {}
        self._param_info[p_name]['type'] = p_type
        self._param_info[p_name]['info'] = None

        if p_type == 'uniform':
            p_info = [float(p_info[0]), float(p_info[1])]
        elif p_type == 'equals':
            p_info = [p_info[0]]
        elif p_type == 'static':
            p_info = [float(p_info[0])]
        else:
            err_msg = 'Parameter type of [{}] is not supported'
            err_msg = err_msg.format(p_type)
            raise ValueError(err_msg)
 
        self._param_info[p_name]['info'] = p_info

    def _set_potential_pair_type(self,params):
        element1  = params[0].strip()
        element2  = params[1].strip()
        pair_type = params[2].strip()
        pair_str  = "{}{}".format(element1,element2)
        self._pair_type[pair_str] = {}
        self._pair_type[pair_str]['type'] = pair_type
        self._pair_type[pair_str]['param'] = {}

    def _set_potential_3body_type(self,params):
        element1 = params[0].strip()
        element2 = params[1].strip()
        element3 = params[2].strip()
        pot_type = params[3].strip()
        str_3body = "{}{}{}".format(element1,element2,element3)
        self._pair_type[str_3body] = {}
        self._pair_type[str_3body]['type'] = pot_type
        self._pair_type[str_3body]['param'] = {}

    def get_param_list(self):
        self.param_list = []
        for element in self.elements:
            self.param_list.append('chrg_{}'.format(element,self.elements[element]['charge']))
        for pair in self.pair:
            for param in self.pair[pair]['param']:
                self.param_list.append('p_{}_{}'.format(pair,param))
        return sorted(self.param_list)

    def get_param_value(self,param_name):
        pot_type = param_name.split('_')[0]
        if pot_type == 'p':
            pot_id   = param_name.split('_')[1]
            pot_param = param_name.split('_')[2]
            return self.pair[pot_id]['param'][pot_param]
        elif pot_type == 'chrg':
            pot_id   = param_name.split('_')[1]
            return self.elements[pot_id]['charge']
          
class QoiConfigFile:
    """Class which processes the pyposmat.qoi file"""
    def __init__(self, 
                 fname_config = "pyposmat.qoi",
                 is_read = True):
        """default constuctor

        Args:
            fname_config (str, optional): filename of the pyposmat.qoi file.
                default is pyposmat.qoi
            is_read (bool, optional): if set to true, it will read the
                configuration file.  default is True.
        """
        self._fname_config = fname_config
        self._qoi_names = []
        self._qoi_info = None
        if is_read == True:
            self.read()

    @property
    def qoi_ref_values(self):
        """:obj:`dict` of :obj:'float': dictionary where the key is the 
        qoi name and the value"""
        qoi_info = {}
        for k,v in self._qoi_info.items():
            qoi_info[k] = v['target']
        return qoi_info

    @property
    def qoi_info(self):
        return self._qoi_info

    @qoi_info.setter
    def qoi_info(self, val):
        self._qoi_info = val

    @property
    def qoi_names(self):
        """:obj:`list` of :obj:`str`: list of qoi names"""
        return self._qoi_names

    def read(self):
        #self.qoi_info = {}
        #self.qoi_info['qoi'] = {}
        f = open(self._fname_config)
        lines = f.readlines()
        f.close()
        for line in lines:
            line = line.strip()
            if line.startswith('#') or line == '':
                pass
            else:
                keyword = line.split('=')[0].strip()
                params  = line.split('=')[1].strip()
                if keyword == 'define_qoi':
                    self._define_qoi(params)
                elif keyword == 'qoi_target':
                    self._define_qoi_targets(params)
                else:
                    msg_err = "unknown qoi_keyword({})"
                    msg_err = msg_err.format(keyword)
                    raise ValueError(msg_err)

    def _define_qoi(self,params):
        qoi_info = params.split(',')
                
        qoi_name       = qoi_info[0]
        qoi_variable   = qoi_info[1].strip()
        qoi_structures = [qoi_info[i].strip() for i in range(2,len(qoi_info))]

        if self._qoi_info is None:
            self._qoi_info = {}
        self._qoi_names.append(qoi_name)
        self._qoi_info[qoi_name] = {}
        self._qoi_info[qoi_name]['variable']  = qoi_variable
        self._qoi_info[qoi_name]['structure'] = qoi_structures                

    def _define_qoi_targets(self,params):
        qoi_list     = params.split(',')
        assert len(qoi_list),4

        qoi_name     = qoi_list[0].strip()
        qoi_target   = float(qoi_list[1].strip())
        qoi_weight   = float(qoi_list[2].strip())
        qoi_norm     = float(qoi_list[3].strip())

        if self._qoi_info is None:
            self._qoi_info = {}
        self._qoi_info[qoi_name]['target'] = qoi_target
        self._qoi_info[qoi_name]['weight'] = qoi_weight
        self._qoi_info[qoi_name]['norm']   = qoi_norm

#TODO: Deprecated remove
def calculate_bulk_modulus(c11, c12):
  K = (c11+2*c12)/3.
  return K

#TODO: Deprecated remove
def calculate_shear_modulus(c11, c12):
  G = (c11-c12)/2.
  return G

class ParameterSampler:
    def __init__(self, param_names, param_dict):
        self._param_names = param_names
        self._param_dict = param_dict

    def generate_parameter_set(self):
        raise NotImplementedError

class KdeSampler(ParameterSampler):
    def __init__(self, param_names, param_dict):
        ParameterSampler.__init__(self,param_names,param_dict)

class UniformSampler(ParameterSampler):
    def __init__(self, param_names, param_dict):
        ParameterSampler.__init__(self,param_names,param_dict)
        
    def generate_parameter_sets(self):
        pass

