# -*- coding: utf-8 -*-
"""This module contains the pyposmat engine for parameterization"""
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2016,2017,2018"
__license__ = "Simplified BSD License"
__version__ = "1.0"

import time,sys
from collections import OrderedDict
import scipy.stats
from pypospack.pyposmat import PyposmatEngine
from pypospack.pyposmat import PyposmatDataFile
from pypospack.task.lammps import LammpsSimulationError
from pypospack.task.task_manager import PypospackTaskManagerError

class PyposmatMonteCarloSampler(PyposmatEngine):
    def __init__(self,
            filename_in='pypospack.config.in',
            filename_out='pypospack.results.out',
            base_directory=None):
        assert isinstance(filename_in,str)
        assert isinstance(filename_out,str)
        assert type(base_directory) in [str,type(None)]

        PyposmatEngine.__init__(self,
                filename_in=filename_in,
                filename_out=filename_out,
                base_directory=base_directory,
                fullauto=False)

        self.pyposmat_data_in_filename = None
        self.pyposmat_data_out_filename = filename_out

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
        
        time_start = time.time()
        _n_errors = 0

        _constrained_parameter_names = []
        for p in self.parameter_names:
            if p not in self.free_parameter_names:
                _constrainted_parameter_names.append(p)

        for i_sample in range(n_samples):
            # generate parameter set
            _parameters = OrderedDict([(p,None) for p in self.parameter_names])
            for p in self.free_parameter_names:
                _parameters[p] = _rv_generators[p].rvs(size=1)[0]

            for p in _constrained_parameter_names:
                _str_eval = str(_param_dist_def[p][1])
                for fp in free_parameter_names:
                    if fp in _str_eval:
                        _str_eval = _str_eval.replace(fp,str(_parameters[fp]))
                _parameters[p] = eval(_str_eval)

            try:
                _results = self.evaluate_parameter_set(parameters=_parameters)
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
                    time_total = time_end-time_start
                    avg_time = n_samples_completed/time_total
                    _str_msg = '{} samples completed in {:.4f}s. Avg_time = {:.4f}. n_errors = {}'.format(
                        n_samples_completed,
                        time_total,
                        avg_time,
                        _n_errors)
                    print(_str_msg)
                    sys.stdout.flush()

    def run_kde_sampling(self,n_samples,filename_in):    
        _rv_generators = OrderedDict()
        for p,definiton in self.free_parameter_names.items():
            distribution_type = definition[0]
            if distribution_type == 'uniform':
                _a = definition[1]['a']
                _b = definition[1]['a']
                _loc = _a
                _scale = _b-_a
                _rv_generators[p] = scipy.stats.uniform(loc=_loc,scale=_scale)
            else:
                raise ValueError('unknown distribution type: {}'.format(
                    distribution_type))

        engine.pyposmat_datafile_out.write_header_section(
                filename=filename_out,
                parameter_names=parameter_names,
                qoi_names=qoi_names,
                error_names=error_names)
        
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

