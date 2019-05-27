import os,shutil
from collections import OrderedDict
from pypospack.task.lammps import LammpsStructuralMinimization

def get_pypospack_root_directory():
    import os
    root_dir = [v.strip() for v in os.environ['PYTHONPATH'].split(':') if v.strip().endswith('pypospack')][0]
    return root_dir

pot_definition = OrderedDict()
pot_definition['potential_type'] = 'eam'
pot_definition['symbols'] = ['Ni']
pot_definition['setfl_filename'] = None
pot_definition['pair_type'] = 'bornmayer'
pot_definition['density_type'] = 'eam_dens_exp'
pot_definition['embedding_type'] = 'eam_embed_eos_rose'
pot_definition['N_r'] = 10000
pot_definition['r_max'] = 10.0
pot_definition['r_cut'] = 10.0
pot_definition['N_rho'] = 10000
pot_definition['rho_max'] = 1000.0
pot_definition['symbols'] = ['Ni']
pot_definition['a0'] = 3.52
pot_definition['lattice_type'] = 'fcc'

a0=3.52
r0=1/(2**0.5)*a0

parameters = OrderedDict()
parameters['p_NiNi_phi0'] = 1.0
parameters['p_NiNi_gamma'] = 2.0
parameters['p_NiNi_r0'] = r0
parameters['d_Ni_rho0'] = 1.0
parameters['d_Ni_beta'] = 4.0
parameters['d_Ni_r0'] = r0
parameters['e_Ni_latticetype'] = 'fcc'
parameters['e_Ni_ecoh'] = -4.45
parameters['e_Ni_B']= 188.
parameters['e_Ni_a0'] = a0

structures = OrderedDict()
structures['structure_directory'] = os.path.join(get_pypospack_root_directory(),'data/Ni_structure_db')
structures['structures'] = OrderedDict()
structures['structures']['Ni_fcc'] = 'Ni_fcc_100_unit.gga.relaxed.vasp'

configuration = OrderedDict()
configuration['task'] = OrderedDict()
configuration['task']['task_name'] = 'Ni_fcc_unit.lmps_min_all'
configuration['task']['task_directory'] = 'Ni_fcc_unit.lmps_min_all'
configuration['task']['task_type'] = 'lmps_min_all'
configuration['potential'] = pot_definition
configuration['parameters'] = parameters
configuration['structure'] = structures


task_name = configuration['task']['task_name']
task_directory = configuration['task']['task_directory']
structure_filename = os.path.join(
        configuration['structure']['structure_directory'],
        configuration['structure']['structures']['Ni_fcc'])
restart=False
fullauto=False

task = LammpsStructuralMinimization(
        task_name = task_name,
        task_directory = task_directory,
        structure_filename = structure_filename)
task.on_init(configuration)
task.on_config(configuration)
task.on_ready(configuration)
task.on_running(configuration)
while task.status is not 'POST':
    task.update_status()
task.on_post(configuration)

print(task.results)
