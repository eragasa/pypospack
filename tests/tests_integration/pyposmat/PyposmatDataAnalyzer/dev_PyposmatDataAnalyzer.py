# -*- coding: utf-8 -*-
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2016,2017"
__license__ = "Simplified BSD License"
__version__ = "1.0"

import os, copy
import numpy as np
from pypospack.pyposmat import PyposmatConfigurationFile
from pypospack.pyposmat import PyposmatDataFile
from pypospack.pyposmat import PyposmatDataAnalyzer
import pypospack.pareto as pareto

if __name__ == "__main__":
    
    import Ni__eam__morse_exp_fs as config
    configuration = PyposmatConfigurationFile()
    configuration.qois = config.qoi_db.qois
    configuration.qoi_constraints = config.qoi_constraints
    configuration.structures = config.structure_db
    configuration.potential = config.potential_formalism
    configuration.sampling_type = config.sampling
    configuration.sampling_distribution = config.parameter_distribution
    configuration.sampling_constraints = config.parameter_constraints
    configuration.write(filename='pypospack.config.in')
    configuration.read(filename='pypospack.config.in')
    
    
    data_directory = 'data__morse_exp_fs'
    pyposmat_data_filename = 'pyposmat.results.4.out'
    pyposmat_configuration_filename = 'pypospack.config.in'  
    data_analyzer = PyposmatDataAnalyzer()
    data_analyzer.read_configuration_file(
            filename=pyposmat_configuration_filename)
    data_analyzer.read_data_file(
            filename=os.path.join(
                data_directory,
                pyposmat_data_filename))
    data_analyzer.filter_performance_requirements()
    data_analyzer.calculate_pareto_set()
    data_analyzer.write_kde_file(
            filename='pyposmat.kde.out')
    exit()
