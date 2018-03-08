import os,shutil
from collections import OrderedDict
from pypospack.task.lammps import LammpsStructuralMinimization
from pypospack.task.lammps import LammpsSimulationError
potential_definition = OrderedDict()
potential_definition['potential_type'] = 'eam'
potential_definition['symbols'] = ['Ni']
potential_definition['setfl_filename'] = "eam_potentials/Ni99.eam.alloy"

structures = OrderedDict()
structures['structure_directory'] = 'structure_db'
structures['structures'] = OrderedDict()
structures['structures']['Ni_fcc'] = 'Ni_fcc_100_unit.gga.relaxed.vasp'
structures['structures']['Ni_hcp'] = 'Ni_hcp_ortho.vasp' 

for sn,sfn in structures['structures'].items():
    configuration = OrderedDict()
    configuration['task'] = OrderedDict()
    configuration['task']['task_name'] = '{}.lmps_min_all'.format(sn)
    configuration['task']['task_directory'] = '{}.lmps_min_all'.format(sn)
    configuration['task']['task_type'] = 'lmps_min_all'
    configuration['potential'] = potential_definition
    configuration['parameters'] = None
    configuration['structure'] = structures

    task_name = configuration['task']['task_name']
    task_directory = configuration['task']['task_directory']
    structure_filename = os.path.join(
        structures['structure_directory'],
        structures['structures'][sn]
    )
    restart=False
    fullauto=False

    lammps_task = LammpsStructuralMinimization(
        task_name = task_name,
        task_directory = task_directory,
        structure_filename = structure_filename)
    lammps_task.on_init(configuration)
    lammps_task.on_config(configuration)
    lammps_task.on_ready(configuration)
    lammps_task.on_running(configuration)
    try:
        while lammps_task.status is not 'POST':
            lammps_task.update_status()
    except LammpsSimulationError as e:
        if str(e) == "Lammps exited with status 1":
           lmps_output_fn = os.path.join(
               task_directory,'lammps.out'
           )
           with open(lmps_output_fn,'r') as f:
               lines = f.readlines()
           print(lines[-1])
    lammps_task.on_post(configuration)

    for k,v in lammps_task.results.items():
        print("{:20} {:+10.6f}".format(k,v))
