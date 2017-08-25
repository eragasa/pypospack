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
    task_obj_config['MgO_NaCl.E_min_all']['class_name'] = 'LammpsSimulation'
    task_obj_config['MgO_NaCl.E_min_all']['chemical_symbols'] = ['Mg','O']
    task_obj_config['MgO_NaCl.E_min_all']['potential_type'] = 'buckingham'
    task_obj_config['MgO_NaCl.E_min_all']['structure_name'] = 'MgO_NaCl'
    task_obj_config['MgO_NaCl.E_min_all']['structure_filename'] = 'MgO_NaCl.vasp'

    task_prereq_config = {}
    task_results = {}
    task_info = {}

    wf = workflow.WorkflowManager(wf_name=wf_name,wf_directory=None)
    wf.initialize_tasks(task_obj_config)
    wf.run()
        

