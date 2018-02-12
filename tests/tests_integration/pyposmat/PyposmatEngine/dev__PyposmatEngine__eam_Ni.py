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
        target=3.508)
Ni_qoi_db.add_qoi(
        qoi_name='Ni_fcc.a0',
        qoi_type='a11_min_all',
        structures=OrderedDict([
            ('ideal','Ni_fcc')
            ]),
        target=3.508)
Ni_qoi_db.add_qoi(
        qoi_name='Ni_fcc.c11',
        qoi_type='c11',
        structures=OrderedDict([
            ('ideal','Ni_fcc')]),
        target=276.)
Ni_qoi_db.add_qoi(
        qoi_name='Ni_fcc.c12',
        qoi_type='c12',
        structures=OrderedDict([
            ('ideal','Ni_fcc')]),
        target=159.)
Ni_qoi_db.add_qoi(
        qoi_name='Ni_fcc.c44',
        qoi_type='c44',
        structures=OrderedDict([
            ('ideal','Ni_fcc')]),
        target=132.)
Ni_structures = OrderedDict()
Ni_structures['structure_directory'] = 'test__PyposmatEngine'
Ni_structures['structures'] = OrderedDict()
Ni_structures['structures']['Ni_fcc'] = 'Ni_fcc.vasp'

Ni_eam_1_potential_definition = OrderedDict()
Ni_eam_1_potential_definition['potential_type'] = 'eam'
Ni_eam_1_potential_definition['symbols'] = ['Ni']
Ni_eam_1_potential_definition['setfl_filename'] = None
Ni_eam_1_potential_definition['pair_type'] = 'morse'
Ni_eam_1_potential_definition['density_type'] = 'eam_dens_exp'
Ni_eam_1_potential_definition['embedding_type'] = 'eam_embed_universal'
Ni_eam_1_potential_definition['N_r'] = 10000
Ni_eam_1_potential_definition['r_max'] = 10.0
Ni_eam_1_potential_definition['r_cut'] = 10.0
Ni_eam_1_potential_definition['N_rho'] = 10000
Ni_eam_1_potential_definition['rho_max'] = 1000.0
Ni_eam_1_potential_definition['symbols'] = ['Ni']
Ni_eam_1_potential_definition['a0'] = 3.52
Ni_eam_1_potential_definition['lattice_type'] = 'fcc'

Ni_eam_parameters = OrderedDict()
Ni_eam_parameters['p_NiNi_D0'] = 0.001114
Ni_eam_parameters['p_NiNi_a'] = 3.429506
Ni_eam_parameters['p_NiNi_r0'] = 2.6813
Ni_eam_parameters['d_Ni_rho0'] = 10.0
Ni_eam_parameters['d_Ni_beta'] = 5.0
Ni_eam_parameters['d_Ni_r0'] = 2.0
Ni_eam_parameters['e_Ni_F0'] = 4.10341782e-3
Ni_eam_parameters['e_Ni_p'] = .80
Ni_eam_parameters['e_Ni_q'] = .20
Ni_eam_parameters['e_Ni_F1'] = -3.09

from pypospack.pyposmat import PyposmatConfigurationFile
Ni_configuration = PyposmatConfigurationFile()
Ni_configuration.qois = Ni_qoi_db.qois
Ni_configuration.potential = Ni_eam_1_potential_definition
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
