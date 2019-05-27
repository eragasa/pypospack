import os
import pypospack.utils
from collections import OrderedDict
from pypospack.pyposmat.engines import PyposmatEngine
from pypospack.pyposmat.data import PyposmatDataFile
# from pypospack.qoi import QoiManager
# from pypospack.task import TaskManager

import yaml
from pypospack.io.filesystem import OrderedDictYAMLLoader

Ni_qoi = OrderedDict()
Ni_qoi['Ni_fcc.E_coh'] = OrderedDict()
Ni_qoi['Ni_fcc.E_coh']['structures'] = {'structure':'Ni_fcc_unit'}
Ni_qoi['Ni_fcc.a1'] = OrderedDict()
Ni_qoi['Ni_fcc.a1']['structures'] = {'structure':'Ni_fcc_unit'}

Ni_potential = OrderedDict()
Ni_potential['potential_type'] = 'eam'
Ni_potential['setfl_filename'] = None
Ni_potential['pair_type'] = 'morse'
Ni_potential['density_type'] = 'eam_dens_exp'
Ni_potential['embedding_type'] = 'eam_embed_universal'
Ni_potential['N_r'] = 10000
Ni_potential['r_max'] = 10.0
Ni_potential['r_cut'] = 10.0
Ni_potential['N_rho'] = 10000
Ni_potential['rho_max'] = 10.0
Ni_potential['symbols'] = ['Ni']

Ni_eam_parameters = OrderedDict()
Ni_eam_parameters['p_NiNi_D0'] = 0.001114
Ni_eam_parameters['p_NiNi_a'] = 3.429506
Ni_eam_parameters['p_NiNi_r0'] = 2.6813
Ni_eam_parameters['d_Ni_rho0'] = 10.0
Ni_eam_parameters['d_Ni_beta'] = 5.0
Ni_eam_parameters['d_Ni_r0'] = 2.0
Ni_eam_parameters['e_Ni_F0'] = 4.10341782e-3
Ni_eam_parameters['e_Ni_p'] = 8.96274624
Ni_eam_parameters['e_Ni_q'] = 8.95940869
Ni_eam_parameters['e_Ni_F1'] = -3.09

_structure_directory = 'test_PypospackEngine'
Ni_structures = OrderedDict()
Ni_structures['Ni_fcc_unit'] = os.path.join(
        _structure_directory,
        'Ni_fcc_unit.gga.relax.vasp')

Ni_pypospack_config = OrderedDict()
Ni_pypospack_config['qoi_info'] = Ni_qoi
Ni_pypospack_config['structures'] = Ni_structures

_pypospack_root_directory = pypospack.utils.get_pypospack_root_directory()
_pyposmat_config_fn = os.path.join(_pypospack_root_directory,'data','MgO_pareto_data','pyposmat.config.in')

engine = PyposmatEngine(
        filename_in = _pyposmat_config_fn,
        filename_out = 'pypospack.config.out')
engine.configure()

