import os,shutil
from collections import OrderedDict

import pypospack.crystal as crystal
import pypospack.io.vasp as vasp
import pypospack.potential as potential

# An example script for running a LAMMPS simulation for the Stillinger-Weber 
# parameterization
# References:
#     http://lammps.sandia.gov/doc/pair_sw.html

Si_sw_potential_definition = OrderedDict()
Si_sw_potential_definition['potential_type'] = 'stillingerweber'
Si_sw_potential_definition['symbols'] = ['Si']

# Potential parameters from Pizzagalli, J Phys Condens Matter 25 (2013) 055801
Si_sw_pizzagalli = OrderedDict()
Si_sw_pizzagalli['SiSiSi_epsilon']=1.04190
Si_sw_pizzagalli['SiSiSi_sigma']=2.128117
Si_sw_pizzagalli['SiSiSi_a']=1.80
Si_sw_pizzagalli['SiSiSi_lambda']=31.0
Si_sw_pizzagalli['SiSiSi_gamma']=1.10
Si_sw_pizzagalli['SiSiSi_costheta0']=-1/3
Si_sw_pizzagalli['SiSiSi_A']=19.0
Si_sw_pizzagalli['SiSiSi_B']=0.65
Si_sw_pizzagalli['SiSiSi_p']=3.5
Si_sw_pizzagalli['SiSiSi_q']=0.5
Si_sw_pizzagalli['SiSiSi_tol']=0

src_structure_dir = 'test_LammpsElasticCalculation'
Si_structure_definition = OrderedDict()
Si_structure_definition['name'] = 'Si_sc'
Si_structure_definition['filename'] = os.path.join(
        src_structure_dir,
        'Si_sc_unit.vasp')

Si_configuration = OrderedDict()
Si_configuration['task'] = OrderedDict()
Si_configuration['task']['task_name'] = 'Si_sc_unit.lmps_min_all'
Si_configuration['task']['task_directory'] = 'Si_sc_unit.lmps_min_all'
Si_configuration['task']['task_type'] = 'lmps_min_all'
Si_configuration['potential'] = Si_sw_potential_definition
Si_configuration['parameters'] = Si_sw_pizzagalli
Si_configuration['structure'] = Si_structure_definition

configuration=Si_configuration

task_name = configuration['task']['task_name']
task_directory = configuration['task']['task_directory']
structure_filename = configuration['structure']['filename']
restart=False
fullauto=False

from pypospack.task.lammps import LammpsElasticCalculation
lammps_task = LammpsElasticCalculation(
        task_name = task_name,
        task_directory = task_directory,
        structure_filename = structure_filename)
lammps_task.on_init(configuration)
lammps_task.on_config(configuration)
lammps_task.on_ready(configuration)
lammps_task.on_running(configuration)
while lammps_task.status != 'POST':
    lammps_task.update_status()
lammps_task.on_post(configuration)

print(lammps_task.results)
