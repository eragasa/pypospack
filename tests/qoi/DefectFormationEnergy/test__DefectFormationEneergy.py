import pytest
import testing_set_Si_sw

from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.task import TaskManager
from pypospack.qoi import QoiDatabase
from pypospack.qoi import QoiManager
from pypospack.qoi import DefectFormationEnergy

def dev____init__():
    testing_set = testing_set_Si_sw.get_testing_set_Si()

    qoi_names = []
    qoi_infos = []
    for qoi_name,qoi_info in testing_set['qoi_db'].qois.items():
        if qoi_info['qoi_type'] == 'E_formation':
            qoi_names.append(qoi_name)
            qoi_infos.append(qoi_info)

    # only do the first one
    o = DefectFormationEnergy(
            qoi_name=qoi_names[0],
            structures=qoi_infos[0]['structures']
    )

def dev__determine_tasks():
    testing_set = testing_set_Si_sw.get_testing_set_Si()

    qoi_names = []
    qoi_infos = []
    for qoi_name,qoi_info in testing_set['qoi_db'].qois.items():
        if qoi_info['qoi_type'] == 'E_formation':
            qoi_names.append(qoi_name)
            qoi_infos.append(qoi_info)

    # only do the first one
    o = DefectFormationEnergy(
            qoi_name=qoi_names[0],
            structures=qoi_infos[0]['structures']
    )
    o.determine_tasks()

    for k,v in o.tasks.items():
        print(k,v)

def dev__defect_calculation():
    print(80*'-')
    print('{:^80}'.format('defect_calculation'))
    print(80*'-')
    simulation_dir = 'rank_test'

    testing_set = testing_set_Si_sw.get_testing_set_Si()

    configuration = PyposmatConfigurationFile()
    configuration.qois = testing_set['qoi_db'].qois
    configuration.structures = testing_set['structure_db']

    qoi_manager = QoiManager(
            qoi_database=configuration.qois,
            fullauto=False
    )
    qoi_manager.configure()
    qoi_manager.determine_tasks()

    print('qoi_manager.tasks')
    print(len('qoi_manager.tasks')*'-')
    for k,v in qoi_manager.tasks.items():
        print(k,v)

    task_manager = TaskManager(base_directory=simulation_dir)
    task_manager.configure(
            tasks = qoi_manager.tasks,
            structures = testing_set['structure_db']
    )
    task_managger.evaluate_tasks()
    
    qoi_manager.calculate_qois(
            task_results=task_manager.results
    )
    qoi_manager.qois

if __name__ == "__main__":
    dev____init__()
    dev__determine_tasks()
    dev__defect_calculation()
