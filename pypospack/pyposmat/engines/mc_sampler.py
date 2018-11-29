# -*- coding: utf-8 -*-
"""This module contains the pyposmat engine for parameterization"""
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2016,2017,2018"
__license__ = "Simplified BSD License"
__version__ = "1.0"

import time,sys,os,copy,shutil,importlib
from collections import OrderedDict
import numpy as np
import scipy.stats
from pypospack.kde import Chiu1999_h
from pypospack.pyposmat.engines import PyposmatEngine
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.task.task_manager import PypospackTaskManagerError
from pypospack.potential import PotentialObjectMap
from pypospack.pyposmat.data import PyposmatLogFile

# necessary exceptions to handle for this class
from numpy.linalg import LinAlgError
from pypospack.exceptions import LammpsSimulationError
from pypospack.exceptions import PyposmatBadParameterError
from pypospack.exceptions import PypospackBadKdeBandwidthType

class PyposmatMonteCarloSampler(PyposmatEngine):
    def __init__(self,
            filename_in='pypospack.config.in',
            # filename_out='pypospack.results.out',
            filename_out='pyposmat.results.out',
            o_log=None,
            mpi_rank=None,
            mpi_size=None,
            base_directory=None):
        """
        Additional attributes are set by the base class pypospack.pyposmat.engines.PyposmatEngine

        Args:
            filename_in (str) - path of the configuration file
            filename_out (str) - path of the output file
            o_log (PyposmatLogFile) - if type(o_log) is a string, then the string is treated as a path in which to log information to.  If type(o_log) is PyposmatLogFile then it is set as an attribute for the refernce.
            mpi_rank (int)
            mpi_size (int)
            base_directory (None)
        Attributes:
            mpi_rank (int) - this is passed in
            mpi_size (int) - this is passed in
            pyposmat_data_in_filename (str) - the path of the datafile to read in
            pyposmat_data_out_filename (str) - the path of the datafile to write simulation results to

        """
        assert isinstance(filename_in,str)
        assert isinstance(filename_out,str)
        assert type(base_directory) in [str,type(None)]

        PyposmatEngine.__init__(self,
                filename_in=filename_in,
                filename_out=filename_out,
                base_directory=base_directory,
                fullauto=False)

        if mpi_rank is None:
            self.mpi_rank = 0
        else:
            self.mpi_rank = mpi_rank

        if mpi_size is None:
            self.mpi_size = 1
        else:
            self.mpi_size = mpi_size

        assert self.mpi_rank < self.mpi_size

        self.mpi_rank=mpi_rank
        self.mpi_size=mpi_size
        self.pyposmat_data_in_filename = None
        self.pyposmat_data_out_filename = filename_out
        self.pyposmat_data_bad_filename = 'pypospack.results.bad'
        
        try:
            self.configure_logger(o_log)
        except TypeError as e:
            m = "Unable to to configure obj_log based on attribute log:{}".format(str(o_log))
            raise TypeError(m)

    def configure_logger(self,o_log=None):
        """
        Configurtion of the log object has different behavior based upon the type passed
        into the argument o_log.  If o_log is PyposmatLogFile, that object will be accessed
        by reference.  A string is assumed to be a filename location.  By default the
        argument for o_log is None, which means logging will go to standard out by means of 
        the print() function.

        Args:
            o_log (str,PyposmatLogFile,None): default: None
        """

        if type(o_log) is PyposmatLogFile:
            self.obj_log = o_log
        elif type(o_log) is str:
            self.obj_log = PyposmatLogFile(filename=o_log)
        elif o_log is None:
            self.obj_log = None
        else:
            m = "log object must be str, PyposmatLogFile, or None"
            raise TypeError(m)
    
    def log(self,str_msg):
        if type(str_msg) is str:
            m = str_msg
        elif type(str_msg) is list:
            m = "\n".join(str_msg)

        if type(self.obj_log) is PyposmatLogFile:
            self.obj_log.write(m)
        print(m)

    def configure_pyposmat_datafile_in(self,filename):
        self.pyposmat_data_in_filename = filename
        self.pyposmat_datafile_in = PyposmatDataFile(filename)

    def configure_pyposmat_datafile_out(self,filename=None):
        if filename is not None:
            assert type(filename) is str
            self.pyposmat_data_out_filename = filename
        self.pyposmat_datafile_out = PyposmatDataFile(filename)

    def read_configuration_file(self,filename=None):
        PyposmatEngine.read_configuration_file(self,filename=filename)
        self.structure_directory = self.configuration.structures['structure_directory']
        self.n_iterations = self.configuration.sampling_type['n_iterations']
        self.parameter_names = [p for p in self.configuration.sampling_distribution]
        self.qoi_names = [k for k in self.configuration.qois]
        self.error_names = ['{}.err'.format(k) for k in self.qoi_names]
        self.parameter_distribution_definition =\
                self.configuration.sampling_distribution
        
        try:
            self.free_parameter_names = [k for k,v in self.parameter_distribution_definition.items() if v[0] != 'equals']
        except KeyError as e:
            print(self.parameter_distribution_definition.items())
            raise
        if self.configuration.sampling_constraints is not None:
            self.parameter_constraints = copy.deepcopy(self.configuration.sampling_constraints)
        else:
            self.parameter_constraints = OrderedDict()

        self.constrained_parameter_names = []
        for p in self.parameter_names:
            if p not in self.free_parameter_names:
                self.constrained_parameter_names.append(p)

    def run_simulations(self,i_iteration,n_samples=None,filename=None):
        i = i_iteration
        _sampling_type = self.configuration.sampling_type[i]['type']
        _n_samples = self.configuration.sampling_type[i]['n_samples']
        
        if self.mpi_rank == 0:
            m = ["R{}: Starting iteration N={}".format(self.mpi_rank,i_iteration)]
            if _sampling_type is "from_file":
                m += ["R{}: Sampling parameters from {}".format(
                    self.mpi_rank,filename)]
            else:
                m += ["R{}: Attemping n_samples={} with sampling_type={}".format(
                    self.mpi_rank,_n_samples,_sampling_type)]
            if filename is not None:
                m += ["R{}: Using file:{}".format(self.mpi_rank,filename)]
            self.log(m)

        if n_samples is not None:
            _n_samples = n_samples

        if _sampling_type == 'parametric':
            self.run_parameteric_sampling(n_samples=_n_samples)
        elif _sampling_type == 'kde':
            if filename is None:
                raise ValueError('cannot do kde sampling with out filename')
            self.run_kde_sampling(n_samples=_n_samples,filename_in=filename)
        elif _sampling_type == 'from_file':
            if filename is None:
                raise ValueError('cannot do filesampling without file')
            self.run_file_sampling(filename)
        else:
            raise ValueError(
                'unknown sampling type:{}'.format(
                    _sampling_type
                )
            )
        if self.mpi_rank == 0:
            print(i_iteration,_n_samples,_sampling_type)

    def write_data_out_header(self):
        self.pyposmat_datafile_out.write_header_section(
                filename=self.pyposmat_data_out_filename,
                parameter_names=self.parameter_names,
                qoi_names=self.qoi_names,
                error_names=self.error_names)

    def run_parameteric_sampling(self,n_samples):

        # create random number generator
        _rv_generators = OrderedDict()
        for p in self.free_parameter_names:
            distribution_type = self.parameter_distribution_definition[p][0]
            if distribution_type == 'uniform':
                _a = self.parameter_distribution_definition[p][1]['a']
                _b = self.parameter_distribution_definition[p][1]['b']
                _loc = _a
                _scale = _b-_a
                _rv_generators[p] = scipy.stats.uniform(loc=_loc,scale=_scale)
            else:
                raise ValueError('unknown distribution type: {}'.format(
                    distribution_type))

        self.write_data_out_header()
        time_start_iteration = time.time()
        _n_errors = 0

        for i_sample in range(n_samples):

            # new OrderedDict to hold in parameter values
            _parameters = OrderedDict([(p,None) for p in self.parameter_names])

            # generate free parameters for ordered dictionary
            for p in self.free_parameter_names:
                _parameters[p] = _rv_generators[p].rvs(size=1)[0]

            # determine parameters determined from equality constraints
            for p in self.constrained_parameter_names:
                _constraint_type = self.parameter_distribution_definition[p][0]
                if _constraint_type == 'equals':

                    # this condition is for fitting EoS for EAM function which
                    # requires a refernce ground state crystal structure
                    if p.endswith('latticetype'):
                        _v = self.parameter_distribution_definition[p][1]
                        _parameters[p] = _v

                    # process evaluation strings
                    elif type(self.parameter_distribution_definition[p][1]) is not list:
                        _str_eval = str(self.parameter_distribution_definition[p][1])

                        # replace string values with numerical values
                        for fp in self.free_parameter_names:
                            if fp in _str_eval:
                                _str_eval = _str_eval.replace(fp,str(_parameters[fp]))

                        # evaluate the string into a float
                        _parameters[p] = eval(_str_eval)
                    else:
                        raise ValueError("oops")

            # additional tasks added here
            for p in self.constrained_parameter_names:
                if self.parameter_distribution_definition[p][0] == 'equals':
                    if type(self.parameter_distribution_definition[p][1]) is list:
                        # required for EAM potentials to calculate dens_max for embedding function
                        if self.parameter_distribution_definition[p][1][0] == 'equilibrium_density':
                            a0 = self.parameter_distribution_definition[p][1][1]
                            latt = self.parameter_distribution_definition[p][1][2]
                            _parameters[p] = self.calculate_equilibrium_density(a0,latt,_parameters)
            
            try:
                # check constraints
                for k,v in self.parameter_constraints.items():
                    _eval_str = v
                    for pn,pv in _parameters.items():
                        _eval_str = _eval_str.replace(pn,str(pv))

                    if eval(_eval_str) is False:
                        m = "failed parameter constraint, {}".format(k)
                        raise PyposmatBadParameterError()
            
                _results = self.evaluate_parameter_set(parameters=_parameters)
            except PyposmatBadParameterError as e:
                _n_errors += 1
            except LammpsSimulationError as e:
                _n_errors += 1
            except PypospackTaskManagerError as e:
                _n_errors += 1
            else:
                self.pyposmat_datafile_out.write_simulation_results(
                        filename=self.pyposmat_data_out_filename,
                        sim_id=i_sample,
                        results=_results)
            finally:
                # print out summaries every 10 solutions
                if (i_sample+1)%10 == 0:
                    n_samples_completed = i_sample+1
                    time_end = time.time()
                    time_total = time_end-time_start_iteration
                    avg_time = time_total/n_samples_completed
                    _str_msg = 'R{}:{} samples completed in {:.4f}s. Avg_time = {:.4f}. n_errors = {}'.format(
                        self.mpi_rank,
                        n_samples_completed,
                        time_total,
                        avg_time,
                        _n_errors)
                    self.log(_str_msg)

    def get_options_kde_bandwidth(self):
        """
        Returns:
            OrderedDict
        """

        kde_options = OrderedDict()
        kde_options['chiu1999'] = OrderedDict()
        kde_options['chiu1999']['reference'] = 'Chiu, S.T. Ann. Stat. 1991, Vol. 19, No 4. 1883-1905'
        kde_options['chiu1999']['doi'] = '10.1214/aos/1176348376'
        kde_options['chiu1999']['description'] = ""
        kde_options['silverman1984'] = OrderedDict()
        kde_options['silverman1984']['reference'] = 'Silverman, B.W. (1986). Density Estimation for Statistics and Data Analysis. London: Chapman & Hall/CRC. p. 48'
        kde_options['silverman1984']['isbn'] = '0-412-24620-1'

    def run_kde_sampling(self,n_samples,filename_in,cluster_id=None,kde_bw_type='chiu1999'):
        """
        cluster_id (int): cluster_id of interest
        kde_bw_type (str): type of bandwidth
        """
        _datafile_in = PyposmatDataFile()
        _datafile_in.read(filename_in)

        if cluster_id is None:
            _free_parameter_names = [str(v) for v in self.free_parameter_names]
            _X = _datafile_in.df[_free_parameter_names].values.T
        else:
            # subselect the dataframe by the cluster_id of interest
            _datafile_in.df = _datafile_in.df.loc[_datafile_in.df['cluster_id'] == cluster_id]
            _X = _datafile_in.df[self.free_parameter_names].loc[_datafile_in.df['cluster_id'] == cluster_id].values.T
            # self.log.write("cluster_id {c} _X.shape={x}".format(c=cluster_id, x=_X.shape))
        #except LinAlgError as e:
        #    raise
        
        # determine bandwidth
        if kde_bw_type == 'Chiu1999':
            _h = Chiu1999_h(_X)
        elif kde_bw_type == 'Silverman1986':
            _h = Silverman1984
        else:
            m = 'kde_bw_type, {}, is not an implemented bandwidth type'
            raise PypospackBadKdeBandwidthType(m)    

        _rv_generator = scipy.stats.gaussian_kde(_X,_h)

        m = 'Chiu1999_h:{}'.format(_h)
        self.log(m)

        self.write_data_out_header()

        time_start_iteration = time.time()
        _n_errors = 0
        
        for i_sample in range(n_samples):
            
            # new OrderedDict to hold in parameter values
            _parameters = OrderedDict([(p,None) for p in self.parameter_names])
            
            # generate free parameters for ordered dictionary
            _free_parameters = _rv_generator.resample(1)
            for i,v in enumerate(self.free_parameter_names):
                _parameters[v] = float(_free_parameters[i,0])

            # determine parameters determined from equality constraints
            for p in self.constrained_parameter_names:
                _constraint_type = self.parameter_disribution_definition[p][0]
                if _constraint_type  == 'equals':
                    
                    # this condition is for fitting EoS for EAM function which
                    # requires a refernce ground state crystal structure
                    if p.endswith('latticetype'):
                        _v = self.parameter_distribution_definition[p][1]
                        _parameters[p] = _v
                    
                    # process evaluation strings
                    elif type(self.parameter_distribution_definition[p][1]) is not list:
                        _str_eval = str(self.parameter_distribution_definition[p][1])

                        # replace string values with numerical values
                        for fp in self.free_parameter_names:
                            if fp in _str_eval:
                                _str_eval = _str_eval.replace(fp,str(_parameters[fp]))

                        # evaluate the string into a float
                        _parameters[p] = eval(_str_eval)
                    else:
                        raise ValueError("oops")

            for p in self.constrained_parameter_names:
                if self.parameter_distribution_definition[p][0] == 'equals':
                    # some EAM potentials have a normalizing equilbirum density 
                    # which have to be determined based upon the parameterization of 
                    # the electron density function
                    if type(self.parameter_distribution_definition[p][1]) is list:
                        if self.parameter_distribution_definition[p][1][0] == 'equilibrium_density':
                            a0 = self.parameter_distribution_definition[p][1][1]
                            latt = self.parameter_distribution_definition[p][1][2]
                            _parameters[p] = self.calculate_equilibrium_density(a0,latt,_parameters)


            try:
                # now we check parameter constraints
                for k,v in self.parameter_constraints.items():
                    _eval_str = v
                    for pn,pv in _parameters.items():
                        _eval_str = _eval_str.replace(pn,str(pv))

                    if eval(_eval_str) is False:
                        s = 'parameter constraint failed, {}'.format(k)
                        raise PyposmatBadParameterError(s)

                _results = self.evaluate_parameter_set(parameters=_parameters)
            except PyposmatBadParameterError as e:
                _n_errors += 1
            except LammpsSimulationError as e:
                _n_errors += 1
            except PypospackTaskManagerError as e:
                _n_errors += 1
            else:
                self.pyposmat_datafile_out.write_simulation_results(
                        filename=self.pyposmat_data_out_filename,
                        sim_id=i_sample,
                        cluster_id=cluster_id,
                        results=_results)
            finally:
                # print out summaries every 10 solutions
                if (i_sample+1)%10 == 0:
                    n_samples_completed = i_sample+1
                    time_end = time.time()
                    time_total = time_end-time_start_iteration
                    avg_time = time_total/n_samples_completed
                    _str_msg = 'R{}:{} samples completed in {:.4f}s. Avg_time = {:.4f}. n_errors = {}'.format(
                        self.mpi_rank,
                        n_samples_completed,
                        time_total,
                        avg_time,
                        _n_errors)
                    self.log(_str_msg)

        d = OrderedDict()
        d['kde_bandwidth'] = OrderedDict()
        d['kde_bandwidth']['type'] = kde_bw_type
        d['kde_bandwidth']['h'] = _h
    
    def run_file_sampling(self,filename_in):

        _datafile_in = PyposmatDataFile(filename=filename_in)
        _datafile_in.read()
        # configure random number generator

        self.write_data_out_header()

        time_start_iteration = time.time()

        _n_errors = 0
        i_sample = 0
        for row in _datafile_in.df.iterrows():
            if self.mpi_rank != i_sample%self.mpi_size:
                i_sample += 1
                continue
            else:
                i_sample += 1
            _parameters = OrderedDict([(p,row[1][p]) for p in self.parameter_names])
            _sim_id = row[1]['sim_id']

            # generate wierd things
            for p in self.constrained_parameter_names:
                if self.parameter_distribution_definition[p][0] == 'equals':
                    if type(self.parameter_distribution_definition[p][1]) is list:
                        if self.parameter_distribution_definition[p][1][0] == 'equilibrium_density':
                            a0 = self.parameter_distribution_definition[p][1][1]
                            latt = self.parameter_distribution_definition[p][1][2]
                            _parameters[p] = self.calculate_equilibrium_density(a0,latt,_parameters)
            try:
                # check constraints
                for k,v in self.parameter_constraints.items():
                    _eval_str = v
                    for pn,pv in _parameters.items():
                        _eval_str = _eval_str.replace(pn,str(pv))
                    if eval(_eval_str) is False:
                        raise PyposmatBadParameterError()

                _results = self.evaluate_parameter_set(parameters=_parameters)
            except PyposmatBadParameterError as e:
                _n_errors += 1
            except LammpsSimulationError as e:
                _n_errors += 1
            except PypospackTaskManagerError as e:
                _n_errors += 1
            else:
                self.pyposmat_datafile_out.write_simulation_results(
                        filename=self.pyposmat_data_out_filename,
                        sim_id=_sim_id,
                        results=_results)
            finally:
                # print out summaries every 10 solutions
                i_sample = i_sample+1
                if (i_sample)%10 == 0:
                    n_samples_completed = i_sample
                    time_end = time.time()
                    time_total = time_end-time_start_iteration
                    avg_time = time_total/n_samples_completed
                    _str_msg = '{} samples completed in {:.4f}s. Avg_time = {:.4f}. n_errors = {}'.format(
                        n_samples_completed,
                        time_total,
                        avg_time,
                        _n_errors)
                    print('rank{}:'.format(self.mpi_rank)+_str_msg)

    def calculate_equilibrium_density(self,a0,latt,parameters):
        _parameters = OrderedDict()
        for k,v in parameters.items():
            if k.startswith('d_'):
                _parameters[k[2:]] = v
            s = k[2:].split('_')[0]
        _potential_type = self.configuration.potential['density_type']
        _symbols = self.configuration.potential['symbols']
        _module_name,_class_name = PotentialObjectMap(
                potential_type=_potential_type)
        try:
            _module = importlib.import_module(_module_name)
            _class = getattr(_module,_class_name)
            _dens_potential = _class(symbols=_symbols)
        except:
            raise

        if latt == 'fcc':
            d = OrderedDict([
                ('1NN',2/(2**0.5)*a0),
                ('2NN',1.000*a0),
                ('3NN',1.225*a0)])
            Z= OrderedDict([
                ('1NN',12),
                ('2NN',6),
                ('3NN',24)])
            rcut = (d['2NN']+d['3NN'])/2.

            rmax = 10.
            r = np.linspace(1,10,5000)*rmax/10
            rho = _dens_potential.evaluate(r,_parameters,rcut)

            rho_e = 0
            for m in Z:
                if d[m] < rcut:
                    rho_e += Z[m]*np.interp(d[m],r,rho[s])

            return rho_e

    def print_structure_database(self):
        m = [
            80*'-',
            '{:^80}'.format('STRUCTURE DATABASE'),
            80*'-',
            'structure_directory:{}'.format(self.structure_directory),
            '',
            '{:^20} {:^20}'.format('name','filename'),
            '{} {}'.format(20*'-',20*'-')
        ]
        m += ['{:20} {:20}'.format(k,v) for k,v in self.structures['structures'].items()]
        self.log(m)

    def print_sampling_configuration(self):
        print(80*'-')
        print('{:^80}'.format('SAMPLING CONFIGURATION'))
        print(80*'-')

        print('{:^10} {:^10} {:^20}'.format(
            'iteration',
            'n_samples',
            'sampling_type'))
        print('{} {} {}'.format(10*'-',10*'-',20*'-'))

        for i in range(self.n_iterations):
            _sample_type = self.configuration.sampling_type[i]['type']
            if _sample_type == 'kde_w_clusters':
                _n_samples = self.configuration.sampling_type[i]['n_samples_per_cluster']
            else:
                _n_samples = self.configuration.sampling_type[i]['n_samples']
            print('{:^10} {:^10} {:^20}'.format(i,_n_samples,_sample_type))

    def print_initial_parameter_distribution(self):
        print(80*'-')
        print('{:80}'.format('INITIAL PARAMETER DISTRIBUTION'))
        print(80*'-')
        for p in self.parameter_distribution_definition:
            if p in self.free_parameter_names:
                str_free = 'free'
                print('{:^20} {:^10} {:^10} {:^10} {:^10}'.format(
                    p,
                    str_free,
                    self.parameter_distribution_definition[p][0],
                    self.parameter_distribution_definition[p][1]['a'],
                    self.parameter_distribution_definition[p][1]['b']))
            else:
                str_free = 'not_free'
                print('{:^20} {:^10}'.format(p,str_free))
