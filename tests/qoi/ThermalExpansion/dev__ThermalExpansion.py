import pytest
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.qoi import QoiDatabase
from pypospack.qoi import ThermalExpansion

MgO_qoi_db = QoiDatabase() 
MgO_qoi_db.add_qoi(
        qoi_name='MgO_NaCl.th_exp',
        qoi_type='thermal_expansion_coefficient',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=0.,
        qoi_options={
            'temperature_min':0,
            'temperature_max':1000,
            'temperature_step':100
            }
        )
print(MgO_qoi_db.to_string())

MgO_structure_db = OrderedDict()
MgO_structure_db['structure_directory'] = '../../../structure_db/MgO_structure_db'
MgO_structure_db['structures'] = OrderedDict()
MgO_structure_db['structures']['MgO_NaCl'] = 'MgO_NaCl_unit.gga.relax.vasp'

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

_qoi_name = 'MgO_NaCl.thermal_expansion'
_structures = OrderedDict([('ideal','MgO_NaCl')])
qoi = ThermalExpansion(
        qoi_name=_qoi_name,
        structures=_structures)
qoi.determine_tasks()

print(80*'-')
print('{:^80}'.format('TASKS DATABASE'))
print(80*'-')

max_len_task_name = 0
max_len_task_type = 0
for k,v in qoi.tasks.items():
    max_len_task_name = max(max_len_task_name,len(k))
    max_len_task_type = max(max_len_task_type,len(v['task_type']))

line_format = \
          "{task_name:^" + str(max_len_task_name) + "} "\
        + "{task_type:^" + str(max_len_task_type) + "} "

for k,v in qoi.tasks.items():
    print(line_format.format(
        task_name = k,
        task_type = v['task_type']))
if __name__ == '__main__':
    pass
