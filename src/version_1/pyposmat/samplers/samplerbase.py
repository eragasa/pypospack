# -*- coding: utf-8 -*-
"""This module contains the pyposmat engine for parameterization"""
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2019"
__license__ = "Simplified BSD License"
__version__ = "1.0"

# 2/17/2019 - EJR
# this object was written to being to abstract the different type of samplers in a
# less obtuse method, consisting of a bunch of if statements which mutated behavior leading to
# a bunch of code which was difficult to maintain.

# some standard imports
import time,sys,os,copy,shutil,importlib
from collections import OrderedDict
import numpy as np
import scipy.stats

# we are wrapping up the pyposmat engine
from pypospack.pyposmat.engines import PyposmatEngine
from pypospack.pyposmat.data import (PyposmatConfigurationFile,
                                     PyposmatDataFile,
                                     PyposmatLogFile)
from pypospack.pyposmat.data import  PyposmatDataFile
from pypospack.pyposmat.data import PyposmatLogFile
from pypospack.potential import PotentialObjectMap

# necessary for kde
from pypospack.kde import Chiu1999_h

# necessary exceptions to handle for this class
from numpy.linalg import LinAlgError
from pypospack.exceptions import LammpsSimulationError
from pypospack.exceptions import PyposmatBadParameterError
from pypospack.exceptions import PypospackBadKdeBandwidthType
from pypospack.exceptions import PypospackTaskManagerError

class PyposmatSampler(PyposmatEngine):
    """ Base Sampling Engine to build other engines upon

    Args:
        config_fn (str): filename of the configuration file
        data_out_fn (str): filename where to output the the simulation results
    Attributes:
        config_fn(str): filename of the configuration file
        data_in_fn(str):filename where to get previous simulation results
        data_out_fn(str):filename where to output the current simulation results
        parameters_fn(str):filename where to output current simulation results
        data_in(:obj:PyposmatDataFile): object for reading in a data file
        data_out(:obj:PyposmatDataFile): object for write out a data file

    """
    def __init__(self,
                 configuration='pyposmat.configuration.yaml',
                 mpi_size=None,
                 base_directory=None):

        assert isinstance(configuration,str) \
            or isinstance(configuration,PyposmatConfiguration)
        assert isinstance(results,str) \
            or isinstance(results,PyposmatData)
        assert
        # check types for the attributes
        assert type(config_fn) is str
        assert type(results_fn) is str
        assert type(data_in_fn) in [type(None),str]
        assert type(o_config) in [type(None),PyposmatConfigurationFile]
        assert type(o_log) in [type(None), PyposmatLogFile]
        assert type(mpi_rank) in [type(None),int]
        assert type(mpi_size) in [type(None),int]
        assert type(base_directory) in [type(None),str]

        super().__init__(
                filename_in=config_fn,
                filename_out=results_fn,
                base_directory=base_directory,
                fullauto=False)

        # default values for mpi_attributes
        self.mpi_rank = 0
        self.mpi_size = 1
        self._configure_mpi_attributes(mpi_rank=mpi_rank,mpi_size=mpi_size)

        # set up necessary filenames
        self.config_fn = config_fn
        self.data_in_fn = None
        self.data_out_fn = results_fn
        self.bad_parameters_fn = bad_parameters_fn

        # data_objects
        self.data_in = None
        self.data_out = None

        # configure log object
        self.obj_log = None
        self._configure_logger(o_log)

        # private attributes
        self._parameter_constraints = None
    def _configure_mpi_attributes(self,mpi_rank,mpi_size):
        # default values, these are set in __init__() but declared here for
        # clarity
        self.mpi_rank = 0
        self.mpi_size = 1

        # we enforce the condition that mpi_rank and mpi_size are both integer types
        # and that the mpi rank_id is less than the total number of mpi_ranks
        if all([type(mpi_rank) is int,type(mpi_size) is int]):
            if mpi_rank < mpi_size:
                self.mpi_rank = mpi_rank
                self.mpi_size = mpi_size

    def _configure_logger(self,o_log=None):
        """configure the logging service

        Configuration of the log object has different behavior based upon the type passed
        into the argument o_log.  If o_log is PyposmatLogFile, that object will be accessed
        by reference.  A string is assumed to be a filename location.  By default the
        argument for o_log is None, which means logging will go to standard out by means of
        the print() function.

        Args:
            o_log (str,PyposmatLogFile,optional): default: None
        Raises:
            TypeError
        """

        assert type(o_log) in [type(None),str,PyposmatLogFile]

        if type(o_log) is PyposmatLogFile:
            self.obj_log = o_log
        elif type(o_log) is str:
            self.obj_log = PyposmatLogFile(filename=o_log)
        elif type(o_log) is type(None):
            self.obj_log = PyposmatLogFile()
        else:
            m = "o_log must be str, PyposmatLogFile, or None"
            raise TypeError(m)

    @property
    def n_iterations(self):
        if type(self.configuration) is not type(None):
            return self.configuration.sampling_type['n_iterations']
        else:
            return None

    @property
    def parameter_names(self):
        if type(self.configuration) is not type(None):
            return self.configuration.parameter_names
        else:
            return None

    @property
    def qoi_names(self):
        if type(self.configuration) is not type(None):
            return self.configuration.qoi_names
        else:
            return None

    @property
    def error_names(self):
        if type(self.configuration) is not type(None):
            return self.configuration.error_names
        else:
            return None

    @property
    def parameter_distribution_definition(self):
        if type(self.configuation) is not type(None):
            return self.configuration.sampling_distribution
        else:
            return None

    @property
    def free_parameter_names(self):
        if type(self.configuration) is not type(None):
            return self.configuration.free_parameter_names
        else:
            return None

    @property
    def parameter_constraints(self):
        if type(self.configuration) is not type(None):
            if type(self._parameter_constraints) is type(None):
                return self.configuration.sampling_constraints
            else:
                return None
        else:
            return None

    @property
    def constrained_parameter_names(self):
        if type(self.configuration) is not type(None):
            return [p for p in self.parameter_names if p not in self.free_parameter_names]
        else:
            return None

    def log(self,str_msg):
        """log message to log file

        Args:
            str_msg (str,list):

        Raises:
            TypeError: If type(str_msg) not either a :obj:str or a :obj:list of :obj:str
        """

        assert type(str_msg) in [str,list]
        if type(str_msg) is list:
            assert all([type(v) is str for v in str_msg])

        self.obj_log.write(m)
        if type(str_msg) is str:
            m = str_msg
        elif type(str_msg) is list:
            m = "\n".join(str_msg)
        else:
            m = "str_msg must be either be a str or a list of str"
            raise TypeError(m)
        self.obj_log.write(m)

    def read_configuration_file(self,filename=None):
        """read the pyposmat configuration file

        This method overrides the inherited method.

        Args:
            filename(str,optional):path of the filename.  If the filename is
        not specified, then the method will run using the class attribute, `config_fn`

        Returns:
            Nothing returned

        Raises:
            TypeError
        """


        # In the previous iteration, this set a bunch of public attributes.  I
        # have reimplemented them as properties because it is much easier for an
        # external developer to understand property implemntation rather than search
        # for a property which maybe mutated.
        # -- EJR, 2/17/2019

        assert type(filename) in [type(None),str]

        if type(filename) is type(None):
            _filename = self.config_fn
        elif type(filename) is str:
            _filename = filename
        else:
            m = "filename must either be a str or NoneType"
            raise(TypeError(m))

        super().read_configuration_file(filename=_filename)

    def configure_pyposmat_datafile_in(self,filename=None):
        """ configures the data_in attribute

        Args:
            filename(str): path of the input file to be used
        """
        assert type(filename) in [type(None),str]

        if type(filename) is str: self.data_in_fn = filename
        _filename = self.data_in_fn
        self.data_in = PyposmatDataFile(filename=_filename)

    def configure_pyposmat_datafile_out(self,filename):
        """ configures the data_out attribute

        Args:
            filename(str): path of the output file to be used
        """

        assert type(filename) in [type(None),str]

        if type(filename) is str: self.data_out_fn = filename
        _filename = self.data_out_fn
        self.data_out = PyposmatDataFile(filename=_filename)

    def initalize_sampler(self):
        raise NotImplementedError

    def generate_free_parameters(self):
        """ stub implementation which needs to be overrided by the inheriting class"""
        free_parameters = OrderedDict()
        for p in self.free_parameter_names:
            free_parameters[p] = 0.
        return free_parameters

    def enforce_parameter_equality_constraints(self,free_parameters):
        constrained_parameters = OrderedDict()
        for p in self.constrained_parameter_names:
            _constraint_type =self.parameter_distribution_definition[p][0]
            if _constraint_type == 'equals':

                if p.endswith('latticetype'):
                    constrainted_parameters[p] = self.parameter_distribution[p][1]

                # evaluate the strings
                elif type(self.parameter_distribution_definition[p][1]) is not list:

                    # get the string to evaluate
                    s = str(self.parameter_distribution_definition[p][1])

                    # replace string values with numerical values
                    for fp in self.free_parameter_names:
                        if fp in s:
                            s = s.replace(fp,str(free_parameters[fp]))

                    # the string can now be evaluated as a float
                    constrainted_parameters[p] = eval(s)

    def enforce_parameter_inequality_constraints(self,parameters):

        # evaluation string
        for k,v in self.parameter_constraints.items():
            eval_str = v
            for pn,pv in parameters.items():
                eval_str = eval_str.replace(pn,str(pv))

            if not eval(eval_str):
                raise PyposmatBadParameterError()

    def run_simulations(self,i_iteration,n_samples=None,filename=None):
        """ base method to override

        """

        assert type(i_iteration) is int
        assert type(n_samples) in [type(None),int]
        assert type(filename) in [type(None),str]


        # define some convenience local variables for readability
        i = i_iteration
        if n_samples is not None:
            _n_samples = self.configuration.sampling_type[i]['n_samples']
        else:
            _n_samples = n_samples

        _sampling_type = self.configuration.sampling_type[i]['type']
        if filename is not None:
            _filename = self.configuration.sampling_type[i][n_samples]
        else:
            pass
