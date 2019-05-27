from collections import OrderedDict
from pypospack.task import TaskManager

# this script is for the development of the TaskManager for the MgO

# the QoiDatabase needs to be done first
from pypospack.qoi import QoiDatabase, QoiManager

import time

class PypospackTaskManagerError(Exception): pass
class TaskManager__TestHarness(TaskManager):
    def evaluate_tasks__testharness(self,parameters,potential):
        _sleep_time = 0.1

        self.results = OrderedDict()
        _configuration = OrderedDict()
        _configuration['potential'] = potential
        _configuration['parameters'] = parameters

        print('evaluating__potential...')
        for k,v in _configuration['potential'].items():
            print(k,v)
        print('evaluating__parameters...')
        for k,v in _configuration['parameters'].items():
            print(k,v)
        def all_simulations_finished(obj_Task):
            _statuses = [o_task.status in ['FINISHED','ERROR']
                    for o_task in obj_Task.values()]
            return all(_statuses)

        _start_time = time.time()
        while not all_simulations_finished(self.obj_Task):
            # iterate over each task, and try to progress the status
            # INIT -> CONFIG
            # CONFIG -> READY
            # READY -> RUNNING
            # RUNNING -> POST
            # POST -> FINISHED
            _max_time_per_simulation = 100
            _time_elapsed = time.time() - _start_time
            if _time_elapsed > _max_time_per_simulation:
                for k_task,o_task in self.obj_Task.items():
                    # kill off process
                    # https://www.programcreek.com/python/example/11892/os.getpgid
                    # https://stackoverflow.com/questions/4789837/how-to-terminate-a-python-subprocess-launched-with-shell-true/4791612#4791612
                    # https://www.codeday.top/2017/06/28/25301.html
                    try:
                        pid = o_task.process.pid
                        pgid = os.getpgid(pid)
                        if pgid == pid:
                            os.killpg(pgid,signal.SIGKILL)
                        else:
                            os.kill(pgid,signal.SIGKILL)
                    except: pass
                raise PypospackTaskManagerError('simulation time exceeded')
            for k_task,o_task in self.obj_Task.items():
                assert isinstance(o_task.configuration,OrderedDict)
                o_task.update_status()
                if o_task.status == 'INIT':
                    _configuration = OrderedDict()
                    _configuration['potential'] = potential
                    _configuration['parameters'] = parameters
                    if 'bulk_structure' in self.tasks[k_task]:
                        _structure_name = self.tasks[k_task]['bulk_structure']
                        _structure_filename = os.path.join(
                                self.structures['structure_directory'],
                                self.structures['structures'][_structure_name])
                        _configuration['bulk_structure'] = _structure_name
                        _configuration['bulk_structure_filename'] = \
                                _structure_filename
                    print(k_task)
                    for k,v in _configuration.items():
                        print("\t",k,v)

                    o_task.on_init(configuration=_configuration)
                elif o_task.status == 'CONFIG':
                    print('{}:CONFIG'.format(k_task))
                    for k,v in self.results.items():
                        print('\t',k,v)
                    o_task.on_config(
                            configuration=_configuration,
                            results=self.results)
                elif o_task.status == 'READY':
                    print('{}:READY'.format(k_task))
                    o_task.on_ready(results=self.results)
                elif o_task.status == 'RUNNING':
                    print('{}:RUNNING'.format(k_task))
                    print('\t',o_task.process)
                    o_task.on_running()
                elif o_task.status == 'POST':
                    print('{}:POST'.format(k_task))
                    o_task.on_post()
                    _results = o_task.results
                    for k,v in o_task.results.items():
                        print('\t',k,v)
                        self.results[k] = v
                elif o_task.status == 'FINISHED':
                    o_task.on_finished()
                elif o_task.status == 'ERROR':
                    raise ValueError
                else:
                    raise ValueError
             
            time.sleep(_sleep_time)

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
_qoidb_OrderedDict['MgO_NaCl.E_fr_a'] = OrderedDict()
_qoidb_OrderedDict['MgO_NaCl.E_fr_a']['qoi_name'] = 'MgO_NaCl.fr_a'
_qoidb_OrderedDict['MgO_NaCl.E_fr_a']['qoi_type'] = 'E_formation_defect'
_qoidb_OrderedDict['MgO_NaCl.E_fr_a']['structures'] = OrderedDict()
_qoidb_OrderedDict['MgO_NaCl.E_fr_a']['structures']['defect'] = 'MgO_NaCl_fr_a'
_qoidb_OrderedDict['MgO_NaCl.E_fr_a']['structures']['ideal'] = 'MgO_NaCl'
_qoidb_OrderedDict['MgO_NaCl.E_fr_a']['target'] = 6.0

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

print(80*'-')
print('{:^80}'.format('QOIS'))
print(80*'-')
for k,v in qoimanager.qois.items():
    print(k,v)

print(80*'-')
print('{:^80}'.format('TASKS'))
for k,v in qoimanager.tasks.items():
    print(k,v)
# structure database
_structuredb_OrderedDict = OrderedDict()
_structuredb_OrderedDict['structure_directory'] = 'structure_db'
_structuredb_OrderedDict['structures'] = OrderedDict()
_structuredb_OrderedDict['structures']['MgO_NaCl'] = 'MgO_NaCl_unit.gga.relax.vasp'
_structuredb_OrderedDict['structures']['MgO_NaCl_fr_a'] = 'MgO_NaCl_333_fr_a.vasp'
_structuredb_OrderedDict['structures']['MgO_NaCl_fr_c'] = 'MgO_NaCl_333_fr_c.vasp'
_structuredb_OrderedDict['structures']['MgO_NaCl_sch'] = 'MgO_NaCl_333_sch.vasp'
_structuredb_OrderedDict['structures']['MgO_NaCl_001s'] = 'MgO_NaCl_001s.vasp'
# potential formalism
import MgO
_potential = MgO.MgO_LewisCatlow['potential']
_parameters = MgO.MgO_LewisCatlow['parameters']
# configure task manager
import os
#_base_directory = 'base_test'
_tasks = qoimanager.tasks
_structures = _structuredb_OrderedDict
task_mgr = TaskManager__TestHarness()
task_mgr.configure(tasks=_tasks,structures=_structures)
task_mgr.evaluate_tasks__testharness(parameters=_parameters,potential=_potential)

print(task_mgr.results)

