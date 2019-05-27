import os
from collections import OrderedDict
from pypospack.qoi import QoiDatabase

import yaml
from pypospack.io.filesystem import OrderedDictYAMLLoader


qoi_db = QoiDatabase()
qoi_db.add_qoi(
        qoi_name='Ni_fcc.E_coh',
        qoi_type='Ecoh_min_all',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=-4.45)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.a0',
        qoi_type='a11_min_all',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=3.52)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.c11',
        qoi_type='c11',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=261.)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.c12',
        qoi_type='c12',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=151.)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.c44',
        qoi_type='c44',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=132.)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.B',
        qoi_type='bulk_modulus',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=188.)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.G',
        qoi_type='shear_modulus',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=101.)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.vac',
        qoi_type='E_formation',
        structures=OrderedDict(
            [
                ('defect','Ni_fcc_vac'),
                ('ideal','Ni_fcc')
            ]
        ),
        target=1.6)
#qoi_db.add_qoi(
#        qoi_name='Ni_fcc.o_int',
#        qoi_type='point_defect',
#        structures=OrderedDict([
#                ('defect','Ni_fcc_o_int'),
#                ('ideal','Ni_fcc')]),
#        target)
#qoi_db.add_qoi(
#        qoi_name='Ni_fcc.t_int',
#        qoi_type='point_defect',
#        structures=OrderedDict([
#                ('defect','Ni_fcc_t_int'),
#                ('ideal','Ni_fcc')]),
#        target)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.100s',
        qoi_type='E_surface',
        structures=OrderedDict(
            [
                ('slab','Ni_fcc_100_s'),
                ('ideal','Ni_fcc_100_unit')
            ]
        ),
        target=1.51e-1)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.110s',
        qoi_type='E_surface',
        structures=OrderedDict(
            [
                ('slab','Ni_fcc_110_s'),
                ('ideal','Ni_fcc_110_unit')
            ]
        ),
        target=1.48e-1)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.111s',
        qoi_type='E_surface',
        structures=OrderedDict(
            [
                ('slab','Ni_fcc_111_s'),
                ('ideal','Ni_fcc_111_unit')
            ]
        ),
        target=1.25e-1)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.esf',
        qoi_type='E_stacking_fault',
        structures=OrderedDict([
                ('defect','Ni_fcc_esf'),
                ('ideal','Ni_fcc_111_unit')]),
        target=7.80e-3)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.isf',
        qoi_type='E_stacking_fault',
        structures=OrderedDict([
                ('defect','Ni_fcc_isf'),
                ('ideal','Ni_fcc_111_unit')]),
        target=1.45e-02)
qoi_db.add_qoi(
        qoi_name='E_Ni_fcc_hcp',
        qoi_type='phase_order',
        structures=OrderedDict([
                ('low','Ni_fcc'),
                ('high','Ni_hcp')]),
        target=0.024)
qoi_db.add_qoi(
        qoi_name='E_Ni_fcc_bcc',
        qoi_type='phase_order',
        structures=OrderedDict([
                ('low','Ni_fcc'),
                ('high','Ni_bcc')]),
        target=0.092)
qoi_db.add_qoi(
        qoi_name='E_Ni_fcc_sc',
        qoi_type='phase_order',
        structures=OrderedDict([
                ('low','Ni_fcc'),
                ('high','Ni_sc')]),
        target=0.600)
qoi_db.add_qoi(
        qoi_name='E_Ni_fcc_dia',
        qoi_type='phase_order',
        structures=OrderedDict([
                ('low','Ni_fcc'),
                ('high','Ni_dia')]),
        target=1.27)

structures = OrderedDict()
structures['structure_directory'] = 'rsrc'
structures['structures'] = OrderedDict()
structures['structures']['Ni_fcc'] = 'Ni_fcc.vasp'

potential_definition = OrderedDict()
potential_definition['potential_type'] = 'eam'
potential_definition['symbols'] = ['Ni']
potential_definition['setfl_filename'] = None
potential_definition['pair_type'] = 'morse'
potential_definition['density_type'] = 'eam_dens_exp'
potential_definition['embedding_type'] = 'eam_embed_universal'
potential_definition['N_r'] = 10000
potential_definition['r_max'] = 10.0
potential_definition['r_cut'] = 10.0
potential_definition['N_rho'] = 10000
potential_definition['rho_max'] = 1000.0
potential_definition['symbols'] = ['Ni']
potential_definition['a0'] = 3.52
potential_definition['lattice_type'] = 'fcc'

potential_parameters = OrderedDict()
potential_parameters['p_NiNi_D0'] = 0.001114
potential_parameters['p_NiNi_a'] = 3.429506
potential_parameters['p_NiNi_r0'] = 2.6813
potential_parameters['d_Ni_rho0'] = 10.0
potential_parameters['d_Ni_beta'] = 5.0
potential_parameters['d_Ni_r0'] = 2.0
potential_parameters['e_Ni_F0'] = 4.10341782e-3
potential_parameters['e_Ni_p'] = .80
potential_parameters['e_Ni_q'] = .20
potential_parameters['e_Ni_F1'] = -3.09

def write_configuration_file(
        filename,
        qois,
        potential_definition,
        structures):
    configuration = PyposmatConfigurationFile()
    configuration.qois = qois
    configuration.potential = potential_definition
    configuration.sampling_distribution = parameter_distribution
    configuration.structures = structures
    configuration.write(filename=filename)

if __name__ == "__main__":
    from pypospack.pyposmat.engines import PyposmatEngine
    from pypospack.pyposmat.data import PyposmatDataFile
    from pypospack.pyposmat.data import PyposmatConfigurationFile
    # from pypospack.qoi import QoiManager
    # from pypospack.task import TaskManager

    filename_config = 'pypospack.config.in'

    write_configuration_file(
            filename=filename_config,
            qois=qoi_db.qois,
            potential_definition=potential_definition,
            structures=structures)

    #Ni_configuration = PyposmatConfigurationFile()
    #Ni_configuration.qois = Ni_qoi_db.qois
    #Ni_configuration.potential = potential_definition
    #Ni_configuration.structures = Ni_structures
    #Ni_configuration.write(filename='pypospack.config.in')
    #Ni_configuration.read(filename='pypospack.config.in')

    engine = PyposmatEngine(
            filename_in = 'pypospack.config.in',
            filename_out = 'pypospack.config.out')
    engine.configure()

    _parameters = potential_parameters
    results = engine.evaluate_parameter_set(parameters=_parameters)

    print(results)
