from pypospack.pyposmat.engines import PyposmatIterativeSampler
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile

import os
import pypospack.utils

_pypospack_root = pypospack.utils.get_pypospack_root_directory()
_simulation_root =  os.path.join(
    _pypospack_root,
    'examples','Ni__eam__born_exp_rose','01_preconditioning_3.5NN')

_pyposmat_configuration_fn = os.path.join(_simulation_root,'data','pyposmat.config.in')
_pyosmat_log_fn = 'pyposmat.log'
_i_iteration = 0
_last_datafile_fn = os.path.join(_simulation_root,'data','pyposmat.kde.1.out')
_new_datafile_fn = 'pyposmat.results.1.out'


engine = PyposmatIterativeSampler(
        configuration_filename=_pyposmat_configuration_fn,
        log_fn='pyposmat.log')
engine.read_configuration_file(
        filename=_pyposmat_configuration_fn)

_data_fn = os.path.join(_simulation_root,'data','pyposmat.results.0.out')
engine.analyze_results(
        i_iteration=0,
        data_fn=_data_fn,
        config_fn=_pyposmat_configuration_fn,
        kde_fn='pyposmat.kde.1.out')
