# -*- coding: utf-8 -*-
"""This module contains the pyposmat engine for parameterization"""
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2017"
__license__ = "Simplified BSD License"
__version__ = "1.0"

import time
import os, shutil, subprocess

import numpy as np
import scipy.stats

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
        fname_config_pyposmat(str)
        fname_config_potential(str)
        fname_config_qoi(str)
        random_seed(int)
        restart(bool)
    """
    def __init__(self,
            fname_config_pyposmat = "pyposmat.config",
            fname_config_potential = "pyposmat.potential",
            fname_config_qoi = "pyposmat.qoi",
            fname_results = "results.out",
            fname_log = "pyposmat.log",
            random_seed = None,restart = False):

        self.supported_qoi = ['a0','a1','a2','a3',
                              'alpha','beta','gamma',
                              'c11','c12','c44',
                              'bulk_modulus',
                              'shear_modulus',
                              'defect_energy',
                              'surface_energy',
                              'stacking_fault_energy',
                              'total_energy']

        self.fname_config_pyposmat = fname_config_pyposmat
        self.fname_config_potential = fname_config_potential
        self.fname_config_qoi = fname_config_qoi
        
        self.restart = restart
        if self.restart is True:
            raise NotImplementedError("Restart method not implemented")

        self._set_random_seed(seed)    

        # determine output
        self._config_results_file(fname_results)
        self._config_log_file(fname_log)

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
        
