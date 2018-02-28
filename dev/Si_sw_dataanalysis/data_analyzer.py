# -*- coding: utf-8 -*-
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2016,2017"
__license__ = "Simplified BSD License"
__version__ = "1.0"

import os, copy
import numpy as np
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatDataAnalyzer
import pypospack.pareto as pareto

if __name__ == "__main__":

    _configuration_filename = 'pyposmat.config.in'
    #--------------------------------------------------------------------------
    # DEFINE CONFIGURATION FILED
    #--------------------------------------------------------------------------
    import Si_sw as config
    configuration = PyposmatConfigurationFile()
    configuration.qois = config.qoi_db.qois
    configuration.qoi_constraints = config.qoi_constraints
    configuration.structures = config.structure_db
    configuration.potential = config.potential_formalism
    configuration.sampling_type = config.sampling
    configuration.sampling_distribution = config.parameter_distribution
    configuration.sampling_constraints = config.parameter_constraints
    #<------------- write configuration file
    configuration.write(filename=_configuration_filename)
    #--------------------------------------------------------------------------
    # TEST CONFIGURATION FILE
    #--------------------------------------------------------------------------
    configuration.read(filename=_configuration_filename)

    #--------------------------------------------------------------------------
    # ANALYZE DATA 
    #--------------------------------------------------------------------------
    data_directory = 'data'
    pyposmat_data_filename = 'pyposmat.results.9.out'
    pyposmat_configuration_filename = 'pyposmat.config.in'
    data_analyzer = PyposmatDataAnalyzer()
    data_analyzer.read_configuration_file(
            filename=pyposmat_configuration_filename)
    data_analyzer.read_data_file(
            filename=os.path.join(
                data_directory,
                pyposmat_data_filename))
    
    #data_analyzer.filter_performance_requirements()
    #data_analyzer.calculate_pareto_set()

    data_analyzer.write_kde_file(
            filename='pyposmat.kde.out')
    exit()
