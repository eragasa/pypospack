# -*- coding: utf-8 -*-
"""This module contains the pyposmat engine for parameterization"""
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2016,2017"
__license__ = "Simplified BSD License"
__version__ = "1.0"

from pypospack.pyposmat import PyposmatEngine

class PyposmatMonteCarloSampler(PyposmatEngine):
    def __init__(self,
            filename_in='pypospack.config.in',
            filename_out='pypospack.results.out',
            base_directory=None):
        assert isinstance(filename_in,str)
        assert isinstance(filename_out,str)
        assert type(base_directory) in [str,type(None)]

        PyposmatEngine.__init__(self,
                filename_in=filename_in,filename_out=filename_out,base_directory=base_directory,fullauto=False)

        #    def __init__(self,
#            filename_in='pypospack.config.in',
#            filename_out='pypospack.results.out',
#            base_directory=None):

#        assert type(filename_in) is str
#        assert type(filename_out) is str
#        assert type(base_directory) in [str,type(None)]
        
#        PyposmatEngine(self,
#            filename_in=filename_in,
#            filename_out=filename_out,
#            base_directory=base_directory)

#    def evaluate_parameter_sets(self):
#        pass
