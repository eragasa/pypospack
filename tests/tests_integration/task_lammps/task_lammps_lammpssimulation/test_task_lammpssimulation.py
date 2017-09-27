import os, shutil, copy
import pypospack.task.lammps as tsk_lammps
import pypospack.potential
import pytest

class TestTaskLammpsElasticCalculation(object):

    @classmethod 
    def setup_class(self):
        self.task_name = 'task_name'
        self.task_directory = 'task_directory'

        self.structure_db = 'rsrc'
        self.structure_filename = 'MgO_NaCl_unit.vasp'

        self.structure_dict = {}
        self.structure_dict['name'] = 'MgO_NaCl_unit'
        self.structure_dict['filename'] = os.path.join(\
                self.structure_db, 
                self.structure_filename)

        self.param_dict = {}
        self.param_dict['chrg_Mg'] = +2.0
        self.param_dict['chrg_O']  = -2.0
        self.param_dict['MgMg_A']   = 0.0 
        self.param_dict['MgMg_rho'] = 0.5
        self.param_dict['MgMg_C']   = 0.0
        self.param_dict['MgO_A']    = 821.6
        self.param_dict['MgO_rho']  = 0.3242
        self.param_dict['MgO_C']    = 0.0
        self.param_dict['OO_A']     = 2274.00 
        self.param_dict['OO_rho']   = 0.1490
        self.param_dict['OO_C']     = 27.88

        self.potential_dict = {}
        self.potential_dict['potential_type'] = 'buckingham'
        self.potential_dict['elements'] = ['Mg','O']
        self.potential_dict['params'] = copy.deepcopy(self.param_dict)

        self.config_dict = {}
        self.config_dict['structure'] = copy.deepcopy(self.structure_dict)
        self.config_dict['potential'] = copy.deepcopy(self.potential_dict)

        if os.path.exists(self.task_directory):
            shutil.rmtree(self.task_directory)

    def class_init(self):
        self.task = tsk_lammps.LammpsElasticCalculation(
                self.task_name, 
                self.task_directory)

    def class_config(self):
        self.task.config(
                structure = self.structure_dict,
                potential = self.potential_dict)

    def class_ready(self):
        self.task.ready()

    def class_run(self):
        self.task.run(is_mpi=False)

    def class_post(self):
        self.task.postprocess()

    def test_init(self):
        self.class_init()
        assert os.path.exists(self.task_directory)
        assert self.task.task_name == self.task_name
        assert self.task.task_directory == os.path.join(\
                os.getcwd(),
                self.task_directory)
        assert self.task.status == 'INIT'

    def test_config(self):
        self.class_init()
        self.class_config()
        assert self.task.status == 'CONFIG'
        assert isinstance(self.task.potential,pypospack.potential.Buckingham)

    def test_ready(self):
        self.class_init()
        self.class_config()
        self.class_ready()
        assert self.task.status == 'READY'

    def test_run(self):
        self.class_init()
        self.class_config()
        self.class_ready()
        self.class_run()
        assert self.task.status == 'POST'

    def test_post(self):
        self.class_init()
        self.class_config()
        self.class_ready()
        self.class_run()
        self.class_post()
if __name__ == "__main__":
    task_name = 'task_name'
    task_directory = 'task_name'

    structure_db = 'rsrc'
    structure_filename = 'MgO_NaCl_unit.vasp'

    structure_dict = {}
    structure_dict['name'] = 'MgO_NaCl_unit'
    structure_dict['filename'] = os.path.join(\
            structure_db, structure_filename)

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

    config_dict = {}
    config_dict['structure'] = {}
    config_dict['structure']['name'] = 'MgO_NaCl_unit'
    config_dict['structure']['filename'] = 'rsrc/MgO_NaCl_unit.vasp'
    config_dict['filename'] = structure_filename
    config_dict['potential'] = {}
    config_dict['potential']['potential_type'] = 'buckingham'
    config_dict['potential']['elements'] = ['Mg','O']
    config_dict['potential']['params'] = param_dict

    print(80*'-')
    print('task.status: null - > INIT')
  
    task = tsk_lammps.LammpsSimulation(\
            task_name = task_name,
            task_directory = task_directory)

    print(80*'-')
    print('task.status: INIT -> CONFIG')
    if task.status == 'INIT':
        task.config(\
                structure=config_dict['structure'],
                potential=config_dict['potential'])

    print(80*'-')
    print('task.status: CONFIG -> READY')
    if task.status == 'CONFIG':
        pre_var_dict = {}
        task.ready(pre_var_dict)

    print(80*'-')
    print('task.status: READY -> RUN')
    if task.status == 'READY': 
        task.run(param_dict,is_mpi=False)
       
    print(80*'-')
    print('task.status: READY -> POST')
    if task.status == 'RUN':
        pass

    print(80*'-')
    print('task.status: POST -> DONE')
    if task.status == 'POST':
        task.post()
        print(task.results)
        print(task.results['tot_energy'])
