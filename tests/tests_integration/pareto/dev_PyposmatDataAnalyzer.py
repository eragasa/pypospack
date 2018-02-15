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
    
    data_directory = 'data'
    pyposmat_data_filename = 'pypospack.results.out'
    pyposmat_configuration_filename = 'pypospack.config.in'  
    data_analyzer = PyposmatDataAnalyzer()
    data_analyzer.read_configuration_file(
            filename=pyposmat_configuration_filename)
    data_analyzer.read_data_file(
            filename=os.path.join(
                data_directory,
                pyposmat_data_filename))
    data_analyzer.calculate_pareto_set()
    data_analyzer.write_kde_file(
            filename=os.path.join(
                data_directory,
                'pyposmat.kde.out'))

    exit()
