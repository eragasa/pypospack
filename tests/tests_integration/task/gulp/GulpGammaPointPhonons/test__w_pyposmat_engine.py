import os
from pypospack.pyposmat.engines import PyposmatEngine
from pypospack.task.gulp import GulpGammaPointPhonons

structure_directory = '../../../../structure_db/MgO_structure_db'

pyposmat_data_dir = 'data'
pyposmat_config_fn = os.path.join(
        pyposmat_data_dir,'pyposmat.config.in')


from collections import OrderedDict
MgO_LewisCatlow_potential_definition = OrderedDict()
MgO_LewisCatlow_potential_definition['potential_type']='buckingham'
MgO_LewisCatlow_potential_definition['symbols'] = ['Mg','O']
MgO_LewisCatlow_potential_definition['parameter_names']=[
    'chrg_Mg','chrg_O',
    'MgMg_A','MgMg_rho','MgMg_C',
    'MgO_A','MgO_rho','MgO_C',
    'OO_A','OO_rho','OO_C']

MgO_LewisCatlow_parameters = OrderedDict()
MgO_LewisCatlow_parameters['chrg_Mg'] = +2.0
MgO_LewisCatlow_parameters['chrg_O']  = -2.0
MgO_LewisCatlow_parameters['MgMg_A']   = 0.0 
MgO_LewisCatlow_parameters['MgMg_rho'] = 0.5
MgO_LewisCatlow_parameters['MgMg_C']   = 0.0
MgO_LewisCatlow_parameters['MgO_A']    = 821.6
MgO_LewisCatlow_parameters['MgO_rho']  = 0.3242
MgO_LewisCatlow_parameters['MgO_C']    = 0.0
MgO_LewisCatlow_parameters['OO_A']     = 2274.00 
MgO_LewisCatlow_parameters['OO_rho']   = 0.1490
MgO_LewisCatlow_parameters['OO_C']     = 27.88

MgO_LewisCatlow = OrderedDict()

MgO_structures = OrderedDict()
MgO_structures['structure_directory'] = structure_directory
MgO_structures['structures'] = OrderedDict()
MgO_structures['structures']['MgO_NaCl_prim'] = 'MgO_NaCl_prim.gga.relax.vasp'

MgO_task_configuration = OrderedDict()
MgO_task_configuration['task'] = OrderedDict()
MgO_task_configuration['task']['task_name'] = 'MgO_NaCl_prim.gulp_phonons_gamma'
MgO_task_configuration['task']['task_directory'] = 'MgO_NaCl_prim.gulp_phonons_gamma'
MgO_task_configuration['task']['task_type'] = 'gulp_phonons_gamma'
MgO_task_configuration['potential'] = MgO_LewisCatlow_potential_definition
MgO_task_configuration['parameters'] = MgO_LewisCatlow_parameters
MgO_task_configuration['structures'] = MgO_structures

def str__results(o_task):
    s = ""

    max_line_len = 80
    max_k_len = max([len(k) for k in o_task.results])
    max_v_len = max_line_len - max_k_len - 1

    str_format = "{:^"+str(max_k_len)+"} {:^"+str(max_v_len)+"}\n"

    for k,v in o_task.results.items():
        s += str_format.format(k,v)

    return s 

if __name__ == '__main__':
    configuration = MgO_task_configuration

    o_task = GulpGammaPointPhonons(
        task_name = MgO_task_configuration['task']['task_name'],
        task_directory = MgO_task_configuration['task']['task_directory'],
        structure_filename = os.path.join(
            configuration['structures']['structure_directory'],
            configuration['structures']['structures']['MgO_NaCl_prim']
            )
        )

    o_task.on_init(configuration)
    o_task.on_config(configuration)
    o_task.on_ready(configuration)
    o_task.on_running(configuration)

    # we need a loot while the task is running, that polls for the task to be done
    while o_task.status is not 'POST':
        o_task.update_status()

    o_task.on_post(configuration)
   
    print(str__results(o_task))


