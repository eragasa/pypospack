import pytest
from collections import OrderedDict
from pypospack.pyposmat import PyposmatConfigurationFile
from pypospack.qoi import QoiDatabase

import MgO

MgO_qoi_db = QoiDatabase() 
MgO_qoi_db.add_qoi(
        qoi_name='MgO_NaCl.fr_a',
        qoi_type='point_defect',
        structures=OrderedDict([
            ('defect','MgO_NaCl_fr_a'),
            ('ideal','MgO_NaCl')]),
        target=10.978)

print(80*'-')
print('{:^80}'.format('QUANTITIES OF INTEREST'))
print(80*'-')
print(MgO_qoi_db.qois)

MgO_structure_db = OrderedDict()
MgO_structure_db['structure_directory'] = 'test__DefectFormationEnergy'
MgO_structure_db['structures'] = OrderedDict()
MgO_structure_db['structures']['MgO_NaCl'] = 'MgO_NaCl_unit.gga.relax.vasp'
MgO_structure_db['structures']['MgO_NaCl_fr_a'] = 'MgO_NaCl_fr_a.vasp'

print(80*'-')
print('{:^80}'.format('STRUCTURE DATABASE'))
print(80*'-')
print('Structure directory:{}'.format(MgO_structure_db['structure_directory']))
print(20*'-' + ' ' + 59*'-')
print('{:^20} {:^59}'.format('structure_name','structure_filename'))
print(20*'-' + ' ' + 59*'-')
for k,v in MgO_structure_db['structures'].items():
    print('{:<20} {:<59}'.format(k,v))

MgO_configuration = PyposmatConfigurationFile()
MgO_configuration.qois = MgO_qoi_db.qois
MgO_configuration.structures = MgO_structure_db

_qoi_name = 'MgO_NaCl_fr_a.defect_energy'
_structures = MgO_configuration.qois['MgO_NaCl.fr_a']['structures']
from pypospack.qoi import DefectFormationEnergy
qoi = DefectFormationEnergy(
        qoi_name=_qoi_name,
        structures=_structures)
qoi.determine_tasks()

print(80*'-')
print('{:^80}'.format('TASKS DATABASE'))
print(80*'-')
for k,v in qoi.tasks.items():
    tsk_name = k
    tsk_type = v['task_type']
    tsk_structure = v['task_structure']
    tsk_requires = v['task_requires']
    
    if v['task_requires'] is None:
        tsk_requires = None
    elif isinstance(v['task_requires'],list):
        tsk_requires = ','.join(v['task_requires'])
    print(tsk_name,tsk_type,tsk_structure,tsk_structure,tsk_requires)

if __name__ == '__main__':
    pass
