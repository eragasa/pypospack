import pytest
import os
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.task import TaskManager
from pypospack.qoi import QoiDatabase,QoiManager
from pypospack.qoi import StackingFaultEnergyCalculation


def str__task_results_to_table(results):
    """
    results: collections.OrderedDict - results from pypospack.task.TaskManager.results
        attribute.
    """
    table_col_names = ['name','value']
    table_rows = []
    for k,v in results.items():
        table_rows.append([str(k),str(v)])

    return '\n'.join([','.join(v) for v in table_rows])

if __name__ == '__main__':
    potential_definition = OrderedDict()
    potential_definition['potential_type'] = 'eam'
    potential_definition['symbols'] = ['Ni']
    potential_definition['setfl_filename'] = os.path.join(
            'potential_db',
            'Mishin-Ni-Al-2009.eam.alloy'
        )
    #potential_definition['setfl_filename'] = os.path.join(
    #        'potential_db',
    #        'Ni99.eam.alloy'
    #    )

    qoi_db = QoiDatabase()
    qoi_db.add_qoi(
            qoi_name='Ni_fcc.esf',
            qoi_type='E_stacking_fault',
            structures=OrderedDict([
                    ('defect','Ni_fcc_esf'),
                    ('ideal','Ni_fcc_111_unit')]),
            target=7.80e-3)
    qoi_db.add_qoi(
            qoi_name='Ni_fcc.isf',
            qoi_type='E_stacking_fault',
            structures=OrderedDict([
                    ('defect','Ni_fcc_isf'),
                    ('ideal','Ni_fcc_111_unit')]),
            target=1.45e-02)
    qoi_db.add_qoi(
            qoi_name='Ni_fcc.usf',
            qoi_type='E_stacking_fault',
            structures=OrderedDict([
                    ('defect','Ni_fcc_usf'),
                    ('ideal','Ni_fcc_111_unit')]),
            target=1.45e-02)

    structure_db = OrderedDict()
    structure_db['structure_directory'] = 'structure_db'
    structure_db['structures'] = OrderedDict()
    structure_db['structures']['Ni_fcc_111_unit'] = 'Ni_fcc_111_unit.gga.relaxed.vasp'
    structure_db['structures']['Ni_fcc_isf'] = 'Ni_fcc_isf.vasp'
    structure_db['structures']['Ni_fcc_esf'] = 'Ni_fcc_esf.vasp'
    structure_db['structures']['Ni_fcc_usf'] = 'Ni_fcc_usf.vasp'

    configuration = PyposmatConfigurationFile()
    configuration.qois = qoi_db.qois
    configuration.structures = structure_db

    print(80*'-')
    print('{:^80}'.format('QUANTITIES OF INTEREST'))
    print(80*'-')
    print(qoi_db.qois)

    print(80*'-')
    print('{:^80}'.format('STRUCTURE DATABASE'))
    print(80*'-')
    print('Structure directory:{}'.format(structure_db['structure_directory']))
    print(20*'-' + ' ' + 59*'-')
    print('{:^20} {:^59}'.format('structure_name','structure_filename'))
    print(20*'-' + ' ' + 59*'-')
    for k,v in structure_db['structures'].items():
        print('{:<20} {:<59}'.format(k,v))


    qoi_manager = QoiManager(
            qoi_database=configuration.qois,
            fullauto=False)
    qoi_manager.configure()
    qoi_manager.determine_tasks()

    print(80*'-')
    print('{:^80}'.format('TASKS DATABASE'))
    print(80*'-')
    for k,v in qoi_manager.tasks.items():
        tsk_name = k
        tsk_type = v['task_type']
        tsk_structure = v['task_structure']
        print(tsk_name,tsk_type,tsk_structure,tsk_structure)

    task_manager = TaskManager(base_directory='rank_test')
    task_manager.configure(
            tasks=qoi_manager.tasks,
            structures=structure_db
        )
    task_manager.evaluate_tasks(
            parameters=potential_definition,
            potential=potential_definition)

    print(80*'-')
    print('{:^80}'.format('TASK RESULTS'))
    print(80*'-')
    print(str__task_results_to_table(task_manager.results))

    print(80*'-')
    print('{:^80}'.format('INTERMEDIATE CALCULATIONS'))
    print(80*'-')

    E_toten_bulk = task_manager.results['Ni_fcc_111_unit.lmps_min_all.toten']
    n_atoms_bulk = task_manager.results['Ni_fcc_111_unit.lmps_min_all.natoms']
    E_coh_bulk = E_toten_bulk/n_atoms_bulk
    E_toten_sf = task_manager.results['Ni_fcc_esf.lmps_min_sf.toten']
    n_atoms_sf = task_manager.results['Ni_fcc_esf.lmps_min_sf.natoms']

    delta_E = E_coh_bulk*n_atoms_sf - E_toten_sf
    print('E_bulk - E_sf={}'.format(delta_E))
    a1 = task_manager.results['Ni_fcc_esf.lmps_min_sf.a11']
    a2 = task_manager.results['Ni_fcc_esf.lmps_min_sf.a22']
    A = a1*a2
    print('A={}'.format(A))
    print('SFE={}'.format(delta_E/A))

    qoi_manager.calculate_qois(
            task_results=task_manager.results)
    qoi_manager.qois
    for k,v in qoi_manager.qois.items():
        print(k,v['qoi_val'],v['qoi_ref'])
