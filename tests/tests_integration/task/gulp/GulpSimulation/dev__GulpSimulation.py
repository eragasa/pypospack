import os,shutil,copy
from collections import OrderedDict
from pypospack.task.gulp import GulpSimulation

MgO_LewisCatlow=OrderedDict()
MgO_LewisCatlow['potential']=OrderedDict()
MgO_LewisCatlow['parameters']=OrderedDict()

MgO_LewisCatlow['potential']['potential_type']='buckingham'
MgO_LewisCatlow['potential']['symbols']=['Mg','O']
MgO_LewisCatlow['potential']['parameters']=[
        'chrg_Mg','chrg_O',
        'MgMg_A','MgMg_rho','MgMg_C',
        'MgO_A','MgO_rho','MgO_C',
        'OO_A','OO_rho','OO_C']
MgO_LewisCatlow['potential']['r_cut'] = 10.0

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

task_configuration = OrderedDict()
task_configuration['task_name'] = 'gulp_task'
task_configuration['task_directory'] = 'gulp_test_directory'
task_configuration['structure_filename'] = os.path.join(
        'test_GulpSimulation','MgO_NaCl_prim.gga.relax.vasp')
task_configuration['restart'] = False

gulp_task = GulpSimulation(
        task_name=task_configuration['task_name'],
        task_directory=task_configuration['task_directory'],
        structure_filename=task_configuration['structure_filename'],
        restart=task_configuration['restart'])

gulp_task.read_structure_file(filename=None)
gulp_task.configure_potential(
        configuration=MgO_LewisCatlow['potential'])
gulp_task.parameters = MgO_LewisCatlow['parameters']
gulp_task.write_gulp_input_file()
gulp_task.run()
