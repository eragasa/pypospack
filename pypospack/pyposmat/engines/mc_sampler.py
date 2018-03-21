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
from pypospack.task.lammps import LammpsSimulationError
from pypospack.task.task_manager import PypospackTaskManagerError
from pypospack.potential import PotentialObjectMap

from numpy.linalg import LinAlgError

class PyposmatBadParameterError(Exception): pass

class PyposmatMonteCarloSampler(PyposmatEngine):
    def __init__(self,
            filename_in='pypospack.config.in',
            filename_out='pypospack.results.out',
            mpi_rank=None,
            mpi_size=None,
            base_directory=None):
        assert isinstance(filename_in,str)
        assert isinstance(filename_out,str)
        assert type(base_directory) in [str,type(None)]

        PyposmatEngine.__init__(self,
                filename_in=filename_in,
                filename_out=filename_out,
                base_directory=base_directory,
                fullauto=False)
        self.mpi_rank=mpi_rank
        self.mpi_size=mpi_size
        self.pyposmat_data_in_filename = None
        self.pyposmat_data_out_filename = filename_out
        self.pyposmat_data_bad_filename = 'pypospack.results.bad'
    def _log(self,str_msg):
        print(str_msg)

    def configure_pyposmat_datafile_in(self,filename):
        self.pyposmat_data_in_filename = filename
        self.pyposmat_datafile_in = PyposmatDataFile(filename)

    def configure_pyposmat_datafile_out(self,filename=None):
        if filename is not None:
            assert type(filename) is str
            self.pyposmat_data_out_filename = filename
        self.pyposmat_datafile_out = PyposmatDataFile(filename)

    def read_configuration_file(self):
        PyposmatEngine.read_configuration_file(self)
        self.structure_directory = self.configuration.structures['structure_directory']
        self.n_iterations = self.configuration.sampling_type['n_iterations']
        self.parameter_names = [p for p in self.configuration.sampling_distribution]
        self.qoi_names = [k for k in self.configuration.qois]
        self.error_names = ['{}.err'.format(k) for k in self.qoi_names]
        self.parameter_distribution_definition =\
                self.configuration.sampling_distribution
        self.free_parameter_names = [k for k,v in self.parameter_distribution_definition.items() if v[0] != 'equals']

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
        if self.mpi_rank == 0:
            print(i_iteration,_n_samples,_sampling_type)

    def run_parameteric_sampling(self,n_samples):
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

        self.pyposmat_datafile_out.write_header_section(
                filename=self.pyposmat_data_out_filename,
                parameter_names=self.parameter_names,
                qoi_names=self.qoi_names,
                error_names=self.error_names)

        time_start_iteration = time.time()
        _n_errors = 0

        for i_sample in range(n_samples):
            # generate parameter set
            _parameters = OrderedDict([(p,None) for p in self.parameter_names])

            # generate free parameters
            for p in self.free_parameter_names:
                _parameters[p] = _rv_generators[p].rvs(size=1)[0]

            # generate constrained parameters
            for p in self.constrained_parameter_names:
                if self.parameter_distribution_definition[p][0] == 'equals':
                    if type(self.parameter_distribution_definition[p][1]) is not list:
                        _str_eval = str(self.parameter_distribution_definition[p][1])
                        for fp in self.free_parameter_names:
                            if fp in _str_eval:
                                _str_eval = _str_eval.replace(fp,str(_parameters[fp]))
                        _parameters[p] = eval(_str_eval)

            # generate wierd things
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

                    try:
                        _is_constraint_ok = eval(_eval_str)
                    except NameError as e:
                        _str = str(e)
                        _regex_str = "name \'d_NiNi_r0\' is not defined"
                        _err_msg = "BadQoiConstraint:\n"
                        _err_msg += "\t{}".format(k)
                        _err_msg += "\t{}".format(v)
                        _err_msg += "\t{}".format(_eval_str)
                        self._log(_err_msg)
                        raise
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
                        sim_id=i_sample,
                        results=_results)
            finally:
                # print out summaries every 10 solutions
                if (i_sample+1)%10 == 0:
                    n_samples_completed = i_sample+1
                    time_end = time.time()
                    time_total = time_end-time_start_iteration
                    avg_time = time_total/n_samples_completed
                    _str_msg = '{} samples completed in {:.4f}s. Avg_time = {:.4f}. n_errors = {}'.format(
                        n_samples_completed,
                        time_total,
                        avg_time,
                        _n_errors)
                    print(_str_msg)
                    sys.stdout.flush()

    def run_kde_sampling(self,n_samples,filename_in):

        # configure random number generator
        _datafile_in = PyposmatDataFile(filename=filename_in)
        _datafile_in.read()
       
        _X = _datafile_in.df[self.free_parameter_names].values.T
        try:
            _h = Chiu1999_h(_X)
        except LinAlgError as e:
            print('filename:{}'.format(filename_in))
            raise
        _rv_generator = scipy.stats.gaussian_kde(_X,_h)
        print('Chiu199_h:{}'.format(_h))
        #_rv_generator = scipy.stats.gaussian_kde(
        #        _datafile_in.df[self.free_parameter_names].values.T)

        self.pyposmat_datafile_out.write_header_section(
                filename=self.pyposmat_data_out_filename,
                parameter_names=self.parameter_names,
                qoi_names=self.qoi_names,
                error_names=self.error_names)

        time_start_iteration_iteration = time.time()
        _n_errors = 0
        for i_sample in range(n_samples):
            # generate parameter set
            _parameters = OrderedDict([(p,None) for p in self.parameter_names])
            _free_parameters = _rv_generator.resample(1)
            for i,v in enumerate(self.free_parameter_names):
                _parameters[v] = float(_free_parameters[i,0])

            # generate constrained parameters
            for p in self.constrained_parameter_names:
                if self.parameter_distribution_definition[p][0] == 'equals':
                    if type(self.parameter_distribution_definition[p][1]) is not list:
                        _str_eval = str(self.parameter_distribution_definition[p][1])
                        for fp in self.free_parameter_names:
                            if fp in _str_eval:
                                _str_eval = _str_eval.replace(fp,str(_parameters[fp]))
                        _parameters[p] = eval(_str_eval)

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
                        sim_id=i_sample,
                        results=_results)
            finally:
                # print out summaries every 10 solutions
                if (i_sample+1)%10 == 0:
                    n_samples_completed = i_sample+1
                    time_end = time.time()
                    time_total = time_end-time_start_iteration
                    avg_time = time_total/n_samples_completed
                    _str_msg = '{} samples completed in {:.4f}s. Avg_time = {:.4f}. n_errors = {}'.format(
                        n_samples_completed,
                        time_total,
                        avg_time,
                        _n_errors)
                    print(_str_msg)

    def run_file_sampling(self,filename_in):

        _datafile_in = PyposmatDataFile(filename=filename_in)
        _datafile_in.read()
        # configure random number generator

        self.pyposmat_datafile_out.write_header_section(
                filename=self.pyposmat_data_out_filename,
                parameter_names=self.parameter_names,
                qoi_names=self.qoi_names,
                error_names=self.error_names)

        time_start_iteration = time.time()
        if self.mpi_rank is None:
            self.mpi_rank = 0
        if self.mpi_size is None:
            self.mpi_size = 1

        _n_errors = 0
        i_sample = 0
        for row in _datafile_in.df.iterrows():
            if self.mpi_rank != row[0]%self.mpi_size:
                continue
            _parameters = OrderedDict([(p,row[1][p]) for p in self.parameter_names])

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
                        sim_id=i_sample,
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
        print(80*'-')
        print('{:^80}'.format('STRUCTURE DATABASE'))
        print(80*'-')
        print('structure_directory:{}'.format(self.structure_directory))
        print('')
        print('{:^20} {:^20}'.format('name','filename'))
        print('{} {}'.format(20*'-',20*'-'))
        for k,v in self.structures['structures'].items():
            print('{:20} {:20}'.format(k,v))

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
            _n_samples = self.configuration.sampling_type[i]['n_samples']
            _sample_type = self.configuration.sampling_type[i]['type']
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
