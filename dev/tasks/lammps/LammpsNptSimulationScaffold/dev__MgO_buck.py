import os,time
from collections import OrderedDict

import pypospack.utils
from pypospack.task  import LammpsNptSimulation
from lammps_npt_scaffold import LammpsNptSimulationScaffold
# definition of the potential
MgO_LewisCatlow = OrderedDict()
MgO_LewisCatlow['potential_definition'] = OrderedDict()
MgO_LewisCatlow['potential_definition']['potential_type'] = 'buckingham'
MgO_LewisCatlow['potential_definition']['symbols'] = ['Mg','O']
MgO_LewisCatlow['parameters'] = OrderedDict()
MgO_LewisCatlow['parameters']['chrg_Mg'] = +2.0
MgO_LewisCatlow['parameters']['chrg_O']  = -2.0
MgO_LewisCatlow['parameters']['MgMg_A']   = 0.0 
MgO_LewisCatlow['parameters']['MgMg_rho'] = 0.5
MgO_LewisCatlow['parameters']['MgMg_C']   = 0.0
MgO_LewisCatlow['parameters']['MgO_A']    = 821.6
MgO_LewisCatlow['parameters']['MgO_rho']  = 0.3242
MgO_LewisCatlow['parameters']['MgO_C']    = 0.0
MgO_LewisCatlow['parameters']['OO_A']     = 2274.00 
MgO_LewisCatlow['parameters']['OO_rho']   = 0.1490
MgO_LewisCatlow['parameters']['OO_C']     = 27.88

# the structures required for the simulation, which would come from the simulation database
MgO_structure_definition = OrderedDict()
MgO_structure_definition['name'] = 'MgO_NaCl_unit'
MgO_structure_definition['filename'] = os.path.join(
        os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'data','MgO_structure_db',
            'MgO_NaCl_unit.vasp')
)



# the information required to configure the task
npt_task_configuration = OrderedDict()
npt_task_configuration['task_name'] = 'MgO_NaCl_unit.npt.T1000'
npt_task_configuration['task_directory'] = 'MgO_NaCl_unit.npt.T1000'
npt_task_configuration['structure_filename'] = MgO_structure_definition['filename']
npt_task_configuration['temperature'] = 1000 # in Kelvin
npt_task_configuration['pressure'] = 0
npt_task_configuration['time_step'] = 0.001 # in picoseconds
npt_task_configuration['time_total'] = 0.001 * 10000
npt_task_configuration['supercell'] = [10,10,10]

if __name__ == "__main__":
    o = LammpsNptSimulationScaffold(
            potential_definition = MgO_LewisCatlow['potential_definition'],
            structures_definition = MgO_structure_definition,
            task_configuration = npt_task_configuration,
            parameters = MgO_LewisCatlow['parameters'])
    o.on_init(o.task_configuration)
    o.on_config(o.task_configuration)

