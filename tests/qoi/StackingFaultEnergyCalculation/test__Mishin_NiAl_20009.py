import pytest
import os
from collections import OrderedDict

import pypospack.utils


from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.task import TaskManager
from pypospack.qoi import QoiDatabase,QoiManager
from pypospack.qoi import StackingFaultEnergyCalculation

NiAl_eam__Mishin2009 = OrderedDict()
NiAl_eam__Mishin2009['potential_type'] = 'eam'
NiAl_eam__Mishin2009['symbols'] = ['Ni,Al']
NiAl_eam__Mishin2009['setfl_filename'] = os.path.join('potential_db','Mishin-Ni-Al-2009.eam.alloy')
NiAl_eam__Mishin2009['reference'] = 'Purja Pun and Y. Mishin (2009).  Phil. Mag., 89(34-36)'

Ni_fcc_isf = OrderedDict()
Ni_fcc_isf['qoi_name'] = 'Ni_fcc.isf'
Ni_fcc_isf['qoi_type'] = 'E_stacking_fault'
Ni_fcc_isf['structures'] = OrderedDict([('defect','Ni_fcc_isf'),('ideal','Ni_fcc_111_unit')])
Ni_fcc_isf['target'] = 100.
Ni_qoi_db = QoiDatabase()
Ni_qoi_db.add_qoi(
        qoi_name=Ni_fcc_isf['qoi_name'],
        qoi_type=Ni_fcc_isf['qoi_type'],
        structures=Ni_fcc_isf['structures'],
        target=Ni_fcc_isf['target'])

Ni_structure_db = OrderedDict()
Ni_structure_db['structure_directory'] = 'structure_db'
Ni_structure_db['structures'] = OrderedDict()
Ni_structure_db['structures']['Ni_fcc_111_unit'] = 'Ni_fcc_111_unit.gga.relaxed.vasp'
Ni_structure_db['structures']['Ni_fcc_isf'] = 'Ni_fcc_isf.vasp'
Ni_structure_db['structures']['Ni_fcc_esf'] = 'Ni_fcc_esf.vasp'
Ni_structure_db['structures']['Ni_fcc_usf'] = 'Ni_fcc_usf.vasp'


test_data = [
    (NiAl_eam__Mishin2009,Ni_qoi_db,Ni_structure_db)
]
test_id =['NiAl_eam__Mishin2009']
@pytest.mark.parametrize("potential,qoi_db,structure_db", test_data,ids=test_id)
def test__structure_files_exist(potential,qoi_db,structure_db):
    structure_dir = structure_db['structure_directory']
    for struture_id,structure_fn in structure_db['structures'].items():
        filename=os.path.join(structure_dir,structure_fn)
        assert os.path.isfile(filename)

@pytest.mark.parametrize("potential,qoi_db,structure_db", test_data,ids=test_id)
def test__calculate_stacking_fault(potential,qoi_db,structure_db):
    simulation_directory = 'rank_test'

    configuration = PyposmatConfigurationFile()
    configuration.qois = qoi_db.qois
    configuration.structures = structure_db

    qoi_manager = QoiManager(qoi_database=configuration.qois,fullauto=False)
    qoi_manager.configure()
    qoi_manager.determine_tasks()

    task_manager = TaskManager(base_directory='rank_test')
    task_manager.configure(tasks=qoi_manager.tasks,structures=structure_db)
    task_manager.evaluate_tasks(parameters=potential,potential=potential)

    qoi_manager.calculate_qois(
            task_results=task_manager.results)
    qoi_manager.qois

def dev__structure_files_exist(structure_db):
    """

    Args:
        structure_db(OrderedDict)
    """
    assert type(structure_db) is OrderedDict

    m = ['determining if structure files exist']
    structure_dir = structure_db['structure_directory']
    
    len_col_1 = max([len(v) for v in ['True','False']])
    len_col_2 = max([len(v) for v in structures_db['structures'].keys()+['structure_id']])
    len_col_3 = max([len(v) for v in structures_db['structures'].values()+['structure_fn']])
 
    row_fmt = []
    for i in [len_col_1,len_col_2,len_col_3]:
        row_fmt += ['{:'+i+'}']
    row_fmt = ' '.join(row_fmt)

    m = [row_fmt.format('exists','structure_id','structure_fn')]
    m += [row_fmt.format(*[v*'=' for v in [len_col_1,len_col_2,len_col_3]])]
    for structure_id,structure_fn in structure_db['structures'].items():
        filename=os.path.join(structure_dir,structure_fn)
        m += [row_fmt.format(os.path.isfile(filename),structure_id,structure_fn)]
    print("\n",join(m))


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
            qoi_name='Ni_fcc.isf',
            qoi_type='E_stacking_fault',
            structures=OrderedDict([
                    ('defect','Ni_fcc_isf'),
                    ('ideal','Ni_fcc_111_unit')]),
            target=1.45e-02)

    structure_db = OrderedDict()
    structure_db['structure_directory'] = 'structure_db'
    structure_db['structures'] = OrderedDict()
    structure_db['structures']['Ni_fcc_111_unit'] = 'Ni_fcc_111_unit.gga.relaxed.vasp'
    structure_db['structures']['Ni_fcc_isf'] = 'Ni_fcc_isf.vasp' 
    configuration = PyposmatConfigurationFile()
    configuration.qois = qoi_db.qois
    configuration.structures = structure_db

    print(80*'-')
    print('{:^80}'.format('QUANTITIES OF INTEREST'))
    print(80*'-')
    for k,v in qoi_db.qois.items():
        print('qoi_name',k)
        print('\tqoi_type:',v['qoi_type'])
        print('\tstructure_defect:',v['structures']['defect'])
        print('\tstructure_ideal:',v['structures']['ideal'])
        print('\ttarget:',v['target'])

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
    E_toten_sf = task_manager.results['Ni_fcc_isf.lmps_min_sf.toten']
    n_atoms_sf = task_manager.results['Ni_fcc_isf.lmps_min_sf.natoms']

    delta_E = E_toten_sf  - E_coh_bulk*n_atoms_sf
    print('E_bulk - E_sf={}'.format(delta_E))
    a1 = task_manager.results['Ni_fcc_isf.lmps_min_sf.a11']
    a2 = task_manager.results['Ni_fcc_isf.lmps_min_sf.a22']
    A = a1*a2
    print('A={}'.format(A))
    print('SFE={}'.format(delta_E/A))

    qoi_manager.calculate_qois(
            task_results=task_manager.results)
    qoi_manager.qois
    for k,v in qoi_manager.qois.items():
        print(k,v['qoi_val'],v['qoi_ref'])
