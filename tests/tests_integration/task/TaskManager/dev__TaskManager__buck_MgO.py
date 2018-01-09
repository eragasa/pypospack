from collections import OrderedDict
from pypospack.task import TaskManager

# this script is for the development of the TaskManager for the MgO

# the QoiDatabase needs to be done first
from pypospack.qoi import QoiDatabase, QoiManager
_qoidb_OrderedDict = OrderedDict()
_qoidb_OrderedDict['MgO_NaCl.Ecoh'] = OrderedDict()
_qoidb_OrderedDict['MgO_NaCl.Ecoh']['qoi_name'] = 'MgO_NaCl.Ecoh'
_qoidb_OrderedDict['MgO_NaCl.Ecoh']['qoi_type'] = 'Ecoh_min_all'
_qoidb_OrderedDict['MgO_NaCl.Ecoh']['structures'] = OrderedDict()
_qoidb_OrderedDict['MgO_NaCl.Ecoh']['structures']['ideal'] = 'MgO_NaCl'
_qoidb_OrderedDict['MgO_NaCl.Ecoh']['target'] = 4.5
_qoidb_OrderedDict['MgO_NaCl.a11'] = OrderedDict()
_qoidb_OrderedDict['MgO_NaCl.a11']['qoi_name'] = 'MgO_NaCl.a11'
_qoidb_OrderedDict['MgO_NaCl.a11']['qoi_type'] = 'a11_min_all'
_qoidb_OrderedDict['MgO_NaCl.a11']['structures'] = OrderedDict()
_qoidb_OrderedDict['MgO_NaCl.a11']['structures']['ideal'] = 'MgO_NaCl'
_qoidb_OrderedDict['MgO_NaCl.a11']['target'] = 5.0
_qoidb_QoiDatabase = QoiDatabase()
for qoi_name,qoi in _qoidb_OrderedDict.items():
    _qoidb_QoiDatabase.add_qoi(
            qoi_name = qoi['qoi_name'],
            qoi_type = qoi['qoi_type'],
            structures = qoi['structures'],
            target = qoi['target'])

qoimanager = QoiManager(qoi_database=_qoidb_QoiDatabase,fullauto=False)
qoimanager.configure()
qoimanager.determine_tasks()

# structure database
_structuredb_OrderedDict = OrderedDict()
_structuredb_OrderedDict['structure_directory'] = 'test__TaskManager'
_structuredb_OrderedDict['structures'] = OrderedDict()
_structuredb_OrderedDict['structures']['MgO_NaCl'] = 'MgO_NaCl_unit.gga.relax.vasp'

# potential formalism
import MgO
_potential = MgO.MgO_LewisCatlow['potential']
_parameters = MgO.MgO_LewisCatlow['parameters']
# configure task manager
import os
_base_directory = 'base_test'
_tasks = qoimanager.tasks
_structures = _structuredb_OrderedDict
task_mgr = TaskManager(base_directory=_base_directory)
task_mgr.configure(tasks=_tasks,structures=_structures)
task_mgr.evaluate_tasks(parameters=_parameters,potential=_potential)
print(task_mgr.results)
