import pytest 
import os
from collections import OrderedDict

import pypospack.utils
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.task import TaskManager
from pypospack.qoi import QoiDatabase
from pypospack.qoi import QoiManager
from pypospack.qoi import PhaseOrderCalculation

potential= OrderedDict()
potential['potential_type'] = 'eam'
potential['symbols'] = ['Ni']
potential['setfl_filename'] = os.path.join(
        pypospack.utils.get_pypospack_root_directory(),
        'data/potentials/Ni__eam/Mishin-Ni-Al-2009.eam.alloy'
)

qoi_db = QoiDatabase()
qoi_db.add_qoi(
        qoi_name='E_Ni_fcc_hcp',
        qoi_type='phase_order',
        structures=OrderedDict([('low','Ni_fcc'),('high','Ni_hcp')]),
        target=0.24
)

structure_db = OrderedDict()
structure_db['structure_directory'] = os.path.join(
        pypospack.utils.get_pypospack_root_directory(),
        'data/Ni_structure_db')
structure_db['structures'] = OrderedDict()
structure_db['structures']['Ni_fcc'] = 'Ni_fcc_100_unit.gga.relaxed.vasp'
structure_db['structures']['Ni_hcp'] = 'Ni_hcp_ortho.vasp'

task_list = OrderedDict()
task_list['Ni_fcc.lmps_min_all'] = OrderedDict([
    ('task_type','lmps_min_all'),('task_structure','Ni_fcc')])
task_list['Ni_hcp.lmps_min_all'] = OrderedDict([
    ('task_type','lmps_min_all'),('task_structure','Ni_hcp')])

def test__structure_files_exist():
    for structure_id,structure_fn in structure_db['structures'].items():
        filename = os.path.join(structure_db['structure_directory'],structure_fn)
        assert os.path.isfile(filename)

class TestQoi(object):
    def __init__(self,qoi_db,structure_db,potential):
       
        self.configuration = None
        self.qoi_manager = None
        self.task_directory = 'rank_test'
        self.task_manager = None

        self.setup_configuration_file(qoi_db=qoi_db,structure_db=structure_db)
        self.check_structure_db(structure_db)
        self.configure_qoi_manager(qoi_database=self.configuration.qois)
        self.check_task_database()
        self.configure_task_manager(task_directory=self.task_directory,
                                    tasks=self.qoi_manager.tasks,
                                    structures=self.configuration.structures)
        self.evaluate_tasks(parameters=potential,
                            potential=potential)
        self.check_task_results(task_results=self.task_manager.results)
        self.calculate_qois()
        self.check_qoi_results()

    def setup_configuration_file(self,qoi_db,structure_db):
        print('setup_configuration_file')
        self.configuration = PyposmatConfigurationFile()
        self.configuration.qois = qoi_db.qois
        self.configuration.structures = structure_db

    def check_structure_db(self,structure_db):
        print('check_structures')
        structure_directory = self.configuration.structures['structure_directory']
        assert os.path.isdir(structure_directory)
        for structure_name,structure_fn in self.configuration.structures['structures'].items():
            print(structure_name,structure_fn)
            assert os.path.isfile(os.path.join(structure_directory,structure_fn))

    def configure_qoi_manager(self,qoi_database):
        self.qoi_manager = QoiManager(qoi_database=qoi_database,fullauto=False)
        self.qoi_manager.configure()
        self.qoi_manager.determine_tasks()

    def check_task_database(self):
        print('{} {} {}'.format('task_name','task_type','task_structure'))
        for k,v in self.qoi_manager.tasks.items():
            task_name = k
            task_type = v['task_type']
            task_structure = v['task_structure']
            print('{} {} {}'.format(task_name,task_type,task_structure))

            assert k in task_list
            assert task_type == task_list[k]['task_type']
            assert task_structure == task_list[k]['task_structure']

        assert len(self.qoi_manager.tasks) == len(task_list)

    def configure_task_manager(self,task_directory,tasks,structures):
        self.task_manager = TaskManager(base_directory=task_directory)
        self.task_manager.configure(tasks=tasks,
                                    structures=structures)
   
    def evaluate_tasks(self,parameters,potential):
        self.task_manager.evaluate_tasks(parameters=parameters,potential=potential)

    def check_task_results(self,task_results):
        table_col_names = ['name','value']
        table_rows = []
        for k,v in task_results.items():
            table_rows.append([str(k),str(v)])

        print('\n'.join([','.join(v) for v in table_rows]))
    
    def calculate_qois(self):
        self.qoi_manager.calculate_qois(task_results=self.task_manager.results)

    def check_qoi_results(self):
        print(self.qoi_manager.qois)
if __name__ == "__main__":

    o = TestQoi(qoi_db=qoi_db,
                structure_db=structure_db,
                potential=potential)
    
