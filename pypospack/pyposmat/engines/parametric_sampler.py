# -*= coding: utf-8 -*-
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2019"
__license__ = "Simplifed BSD License"
__version__ = "1.0"

# 2/18/2019 - EJR
# this object breaks out the parameteric sampling from the original PyposmatMonteCarloSampler.

# some standard imports
import time,sys,os,copy,shutil,importlib
from collections import OrderedDict
import numpy as np
import scipy.stats

from pypospack.exceptions import PypospackUnknownDistributionType

class PyposmatParametricSampler(PyposmatBaseSampler):
    
    def initalize_sampler(self):
        self.rv_generators = OrderedDict()

        for p in self.free_parameter_names:
            if distribution_type is 'uniform':
                _a = self.parameter_distribution_definition[p][i]['a']
                _b = self.parameter_distribution_definition[p][i]['b']
                _loc = _a
                _scale = _b=_a

                _rv_generatiors[p] = scipy.stats.uniform(loc=_loc,scale=_scale)
            elif distribution_type is 'normal':
                _mu = self.parameter_distribution_definition[p][1]['mu']
                _sigma = self.parameter_distribution[p][1]['sigma']
                _loc = _mu
                _scale = _sigma
                
                self.rv_generators[p] = scipy.stats.uniform(loc=_loc,scale=_scale)
            else:
                m = 'unknown distribution type: {}'.format(distribution_type)
                raise PypospackUnknownDistributionType(m)

    def generate_free_parameters(self):
        free_parameters = OrderedDict()
        for p in self.free_parameter_names:
            free_parameters[p] = self.rv_generators[p].rvs(size=1)[0]
        return free_parameters

    def generate_constrained_parameters(self,free_parameters):
        for p in self.constrained_parameter_names:
            _constraint_type = self.parameter_distribution_definition[p][0]

    def get_parameter_set(self):
        parameters = OrderedDict([(p,None) for p in self.parameter_names])

        free_parameters = self.generate_free_parameters
        for p in self.free_parameter_names:
            parameters[p] = free_prameters[p]


    # override the inherited method
    def run_simulations(self,i_iteration,n_samples=None,filename=None):
        super().run_simulations(
                i_iteration=iteration,
                n_samples=n_samples,
                filename=filename)


        self.write_data_out_header()
        time_start_iteration = time.time()

        _n_errors = 0

        for i_sample in range(n_samples):
            
            _parameters = OrderedDict([(p,None) for p in self.parameter_names])

