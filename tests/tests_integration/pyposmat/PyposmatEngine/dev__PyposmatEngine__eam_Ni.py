import os
from collections import OrderedDict
from pypospack.pyposmat import PyposmatDataFile, PyposmatEngine
from pypospack.pyposmat import QoiDatabase
# from pypospack.qoi import QoiManager
# from pypospack.task import TaskManager

import yaml
from pypospack.io.filesystem import OrderedDictYAMLLoader

Ni_qoi_db = QoiDatabase()
Ni_qoi_db.add_qoi(
        qoi_name='Ni_fcc.E_coh',
        qoi_type='Ecoh_min_all',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=4.5)
Ni_qoi_db.add_qoi(
        qoi_name='Ni_fcc.a0',
        qoi_type='a11_min_all',
        structures=OrderedDict([
            ('ideal','Ni_fcc')
            ]),
        target=4.5)

Ni_structures = OrderedDict()
Ni_structures['structure_directory'] = 'test__PyposmatEngine'
Ni_structures['structures'] = OrderedDict()
Ni_structures['structures']['Ni_fcc'] = 'Ni_fcc.vasp'

Ni_eam_1_potential= OrderedDict()
Ni_eam_1_potential['potential_type'] = 'eam'
Ni_eam_1_potential['setfl_filename'] = None
Ni_eam_1_potential['pair_type'] = 'morse'
Ni_eam_1_potential['density_type'] = 'eam_dens_exp'
Ni_eam_1_potential['embedding_type'] = 'eam_embed_universal'
Ni_eam_1_potential['N_r'] = 10000
Ni_eam_1_potential['r_max'] = 10.0
Ni_eam_1_potential['r_cut'] = 10.0
Ni_eam_1_potential['N_rho'] = 10000
Ni_eam_1_potential['rho_max'] = 10.0
Ni_eam_1_potential['symbols'] = ['Ni']

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

from pypospack.pyposmat import PyposmatConfigurationFile
Ni_configuration = PyposmatConfigurationFile()
Ni_configuration.qois = Ni_qoi_db.qois
Ni_configuration.potential = Ni_eam_1_potential
Ni_configuration.structures = Ni_structures
Ni_configuration.write(filename='pypospack.config.in')
Ni_configuration.read(filename='pypospack.config.in')

engine = PyposmatEngine(
        filename_in = 'pypospack.config.in',
        filename_out = 'pypospack.config.out')
engine.configure()

_parameters = Ni_eam_parameters
results = engine.evaluate_parameter_set(parameters=_parameters)

print(results)
