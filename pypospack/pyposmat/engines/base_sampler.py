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
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatLogFile
from pypospack.potential import PotentialObjectMap

# necessary for kde
from pypospack.kde import Chiu1999_h

# necessary exceptions to handle for this class
from numpy.linalg import LinAlgError
from pypospack.exceptions import LammpsSimuationError
from pypospack.exceptions import PyposmatBadParameterError
from pypospack.exceptions import PypospackBadKdeBadwidthType
from pypospack.exceptions import PypospackTaskManagerError

class PyposmatBaseSampler(PyposmatEngine):
    def __init__(self,
            config_fn='pyposmat.config.in',
            data_out_fn='pyposmat.results.out',
            data_in_fn=None
            o_config = None
            o_log = None
            mpi_rank=None
            mpi_size=None,
            base_directory=None):
        
        # check types for the attributes
        assert type(config_fn) is str
        assert type(data_out_fn) is str
        assert (data_in_fn is None) or (type(data_in_fn is str))
        assert (o_config is None) or (type(o_config) is PyposmatConfigurationFile)
        assert (o_log is None) or (type(o_log) is PyposmatLogFile)
        assert (mpi_rank is None) or (type(mpi_rank) is int)
        assert (mpi_size is None) or (type(mpi_size) is int)
        assert (base_directory is None) or (type(base_directory) is str)

        super().__init__(self,
                filename_in=config_fn,
                filename_out=results_fn,
                base_directory=base_directory,
                fullauto=False)

        # default values for mpi_attributes
        self.mpi_rank = 0
        self.mpi_size = 1
        self._configure_mpi_attributes(mpi_rank=mpi_rank,mpi_size=mpi_size)

        # set up necessary filenames
        self.pyposmat_config_in_fn = config_fn
        self.pyposmat_data_in_fn = None
        self.pyposmat_data_out_fn = filename_out
        self.pyposmat_data_bad_parameters = 'pyposmat.bad_parameters.out'

        # configure log object
        self.obj_log = None
        self._configure_logger(o_log)

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

        assert o_log in [type(None),str,PyposmatLogFile)]

        if type(o_log) is PyposmatLogFile: 
            self.obj_log = o_log
        elif type(o_log) is str:
            self.obj_log = PyposmatLogFile(filename=o_log)
        elif type(o_log) is type(None):
            self.obj_log = None
        else:
            m = "o_log must be str, PyposmatLogFile, or None"
            raise TypeError(m)

    def log(self,str_msg):
        """log message to log file

        Args:
            str_msg (str,list):
            
        Raises:
            TypeError: If type(str_msg) not either a :obj:str or a :obj:list of :obj:str
        """

        assert type(str_msg) in [str,list]
        if type(str_msg) is list: assert all([type(v] is str for v in str_msg])

        self.obj_log.write(m)
        if type(str_msg) is str:
            m = str_msg
        elif type(str_msg) is list:
            m = "\n".join(str_msg)
        else:
            m = "str_msg must be either be a str or a list of str"
            raise TypeError(m)
        self.obj_log.write(m)

    def configure_pyposmat_datafile_in(self,filename):

        assert type(filename) is str

        self.pyposmat_data_in = filename
        self.data_in = PyposmatDatafile(filename=filename)
