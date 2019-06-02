"""
example on how to GMM

adapted from
Jake VanderPlass, Python Data Science Handbook
https://jakevdp.github.io/PythonDataScienceHandbook/05.12-gaussian-mixtures.html
"""
import os
import numpy as np
import pandas as pd

# graphics imports
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

# pypospack imports
import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile


    
# pypospack root directory
pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
print(pypospack_root_dir)
config_fn = os.path.join(
        pypospack_root_dir,
        'data','Si__sw__data','pareto_optimization_unconstrained',
        'pyposmat.config.in')
print(config_fn)
data_fn = os.path.join(
        pypospack_root_dir,
        'data','Si__sw__data','pareto_optimization_unconstrained',
        'pyposmat.kde.20.out')
ref_config_fn = os.path.join(
        pypospack_root_dir,
        'data','Si__sw__data','reference_potentials',
        'pyposmat.config.in')
ref_data_fn = os.path.join(
        pypospack_root_dir,
        'data','Si__sw__data','reference_potentials',
        'pyposmat.kde.1.out')

o_config = PyposmatConfigurationFile()
o_config.read(filename=config_fn)

o_data = PyposmatDataFile()
o_data.read(filename=data_fn)
o_data.create_normalized_errors(
        normalize_type='by_qoi_target',
        qoi_targets=o_config.qoi_targets)
print(o_data.df.columns)
