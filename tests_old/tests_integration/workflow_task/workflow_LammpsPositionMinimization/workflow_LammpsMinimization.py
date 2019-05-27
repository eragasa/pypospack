import os,sys
import importlib, copy
import pypospack.workflow as workflow
import pypospack.task as task 
import pypospack.task.lammps as tsk_lammps

import unittest

if __name__ == "__main__":
    wf_name = 'wf_test'
    wf_directory = 'wf_test'

    task_obj_config = {}
    task_obj_config['MgO_NaCl.E_min_all'] = {}
    task_obj_config['MgO_NaCl.E_min_all']['module_name'] = 'pypospack.task.lammps'
    task_obj_config['MgO_NaCl.E_min_all']['class_name'] = 'LammpsMinimizeStructure'
    task_obj_config['MgO_NaCl.E_min_all']['chemical_symbols'] = ['Mg','O']
    task_obj_config['MgO_NaCl.E_min_all']['potential_type'] = 'buckingham'
    task_obj_config['MgO_NaCl.E_min_all']['structure_name'] = 'MgO_NaCl'
    task_obj_config['MgO_NaCl.E_min_all']['structure_filename'] = 'rsrc/MgO_NaCl_unit.vasp'

    param_dict = {}
    param_dict['chrg_Mg'] = +2.0
    param_dict['chrg_O']  = -2.0
    param_dict['MgMg_A']   = 0.0 
    param_dict['MgMg_rho'] = 0.5
    param_dict['MgMg_C']   = 0.0
    param_dict['MgO_A']    = 821.6
    param_dict['MgO_rho']  = 0.3242
    param_dict['MgO_C']    = 0.0
    param_dict['OO_A']     = 2274.00 
    param_dict['OO_rho']   = 0.1490
    param_dict['OO_C']     = 27.88

    task_prereq_config = {}
    task_results = {}
    task_info = {}

    wf = workflow.WorkflowManager(wf_name=wf_name,wf_directory=None)
    wf.param_dict = copy.deepcopy(param_dict)
    wf.initialize_tasks(task_obj_config)
    wf.run()
        

