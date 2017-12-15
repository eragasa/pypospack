import pytest
import os,shutil
from collections import OrderedDict
import pypospack.crystal as crystal
import pypospack.io.vasp as vasp
from pypospack.potential import Potential

from MgO_buck import MgO_LewisCatlow
from MgO_buck import MgO_structures

MgO_LC_configuration = OrderedDict()
MgO_LC_configuration['task'] = OrderedDict()
MgO_LC_configuration['task']['task_name'] = 'MgO_NaCl.lmps_min_all'
MgO_LC_configuration['task']['task_directory'] = 'MgO_NaCl.lmps_min_all'
MgO_LC_configuration['task']['structure_filename'] = os.path.join(
        MgO_structures['structure_db_dir'],
        MgO_structures['MgO_NaCl_unit']['filename'])

MgO_LC_configuration['potential'] = MgO_LewisCatlow['potential']
MgO_LC_configuration['parameters'] = MgO_LewisCatlow['parameters']

MgO_LC_configuration['expected'] = OrderedDict()
MgO_LC_configuration['expected']['task_type'] = 'lmps_min_all'

configuration=MgO_LC_configuration

def cleanup(task_directory):
    if os.path.isdir(task_directory):
        shutil.rmtree(task_directory)

if __name__ == "__main__":
    task_name = configuration['task']['task_name']
    task_directory = configuration['task']['task_directory']
    structure_filename = configuration['task']['structure_filename']
    restart=False
    fullauto=False

    cleanup(task_directory)
    #<--- code setup
    import time
    from pypospack.task.lammps import LammpsStructuralMinimization
    lammps_task = LammpsStructuralMinimization(
            task_name = task_name,
            task_directory = task_directory,
            structure_filename = structure_filename)
    lammps_task.on_init(configuration)
    lammps_task.on_config(configuration)
    lammps_task.on_ready(configuration)
    
    while lammps_task.status != 'POST':
        lammps_task.update_status()
        print(
                lammps_task.status,
                lammps_task.process is None,
                lammps_task.conditions_POST,
                lammps_task.process.poll()
            )
        time.sleep(0.1)

    lammps_task.on_post(configuration)
    print(lammps_task.results)
