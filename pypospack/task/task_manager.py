import os,copy,importlib,time
from collections import OrderedDict

TaskToClassMap = OrderedDict()
TaskToClassMap['lmps_min_all'] = OrderedDict()
TaskToClassMap['lmps_min_all']['module'] = 'pypospack.task.lammps'
TaskToClassMap['lmps_min_all']['class'] = 'LammpsStructuralMinimization' 
TaskToClassMap['lmps_min_pos'] = OrderedDict()
TaskToClassMap['lmps_min_pos']['module'] = 'pypospack.task.lammps'
TaskToClassMap['lmps_min_pos']['class'] = 'LammpsPositionMinimization'
TaskToClassMap['lmps_min_none'] = OrderedDict()
TaskToClassMap['lmps_min_none']['module'] = 'pypospack.task.lammps'
TaskToClassMap['lmps_min_none']['class'] = 'LammpsSinglePointCalculation'
TaskToClassMap['lmps_elastic'] = OrderedDict()
TaskToClassMap['lmps_elastic']['module'] = 'pypospack.task.lammps'
TaskToClassMap['lmps_elastic']['class'] = 'LammpsElasticCalculation'
TaskToClassMap['lmps_npt'] = OrderedDict()
TaskToClassMap['lmps_npt']['module'] = 'pypospack.task.lammps'
TaskToClassMap['lmps_npt']['class']= 'LammpsNptSimulation'
TaskToClassMap['lmps_neb'] = OrderedDict()
TaskToClassMap['lmps_neb']['module'] = 'pypospack.task.lammps'
TaskToClassMap['lmps_neb']['class'] = 'LammpsNebCalculation'

class TaskManager(object):
    """
    Args:
        base_directory(str)
        tasks(collections.OrderedDict)
        structures(collections.OrderedDict)
        fullauto(bool)
    """

    def __init__(self,base_directory=None,tasks=None,structures=None,
            fullauto=True):
        assert any([
                isinstance(base_directory,str),
                type(base_directory in [type(None)])
            ])
        assert any([\
                isinstance(tasks,OrderedDict),
                type(tasks) in [type(None)]
            ])
        assert any([\
                isinstance(structures,OrderedDict),
                type(structures) in [type(None)]
            ])

        self.tasks = None
        self.structures = None
        self.obj_Task = None
        self.base_directory = None
        self.results = None
        if base_directory is None:
            self.base_directory = os.getcwd()
        elif isinstance(base_directory,str):
            self.base_directory = base_directory
            if not os.path.exists(self.base_directory):
                os.mkdir(self.base_directory)
        if structures is not None:
            self.structures = copy.deepcopy(structures)
        if tasks is not None:
            self.tasks = copy.deepcopy(structures)

        fullauto_conditions = [
                fullauto,
                tasks is not None,
                structures is not None,
            ]

        if all(fullauto_conditions):
            self.process_tasks()

    def configure(self,tasks,structures):
        assert isinstance(tasks,OrderedDict)
        assert isinstance(structures,OrderedDict)
        self.tasks = copy.deepcopy(tasks)
        self.structures = copy.deepcopy(structures)
        self.process_tasks()
    
    def evaluate_tasks(self,parameters,potential):
        _sleep_time = 0.1

        self.results = OrderedDict()
        _configuration = OrderedDict()
        _configuration['potential'] = potential
        _configuration['parameters'] = parameters
        def all_simulations_finished(obj_Task):
            _statuses = [o_task.status in ['FINISHED','ERROR']
                    for o_task in obj_Task.values()]
            return all(_statuses)

        while not all_simulations_finished(self.obj_Task):
            # iterate over each task, and try to progress the status
            # INIT -> CONFIG
            # CONFIG -> READY
            # READY -> RUNNING
            # RUNNING -> POST
            # POST -> FINISHED

            for k_task,o_task in self.obj_Task.items():
                assert isinstance(o_task.configuration,OrderedDict)
                o_task.update_status()
                if o_task.status == 'INIT':
                    _config = copy.deepcopy(_configuration)
                    if 'bulk_structure' in self.tasks[k_task]:
                        _config['bulk_structure'] = self.tasks[k_task]['bulk_structure']
                    o_task.on_init(configuration=_config)
                elif o_task.status == 'CONFIG':
                    o_task.on_config()
                elif o_task.status == 'READY':
                    o_task.on_ready(results=self.results)
                elif o_task.status == 'RUNNING':
                    o_task.on_running()
                elif o_task.status == 'POST':
                    o_task.on_post()
                    _results = o_task.results
                    for k,v in o_task.results.items():
                        self.results[k] = v
                elif o_task.status == 'FINISHED':
                    o_task.on_finished()
                elif o_task.status == 'ERROR':
                    raise ValueError
                else:
                    raise ValueError
             
            time.sleep(_sleep_time)

    def process_tasks(self,tasks=None):
        assert any([
                isinstance(tasks,OrderedDict),
                type(tasks)==type(None)
            ])
        
        if isinstance(tasks,OrderedDict):
            self.tasks = copy.deepcopy(tasks)

        if self.obj_Task is None:
            self.obj_Task = OrderedDict()

        for k_task,v_task in self.tasks.items():
            _task_name = k_task
            _task_type = v_task['task_type']
            _task_structure = v_task['task_structure']
            self.add_task(
                    task_name=_task_name,
                    task_type=_task_type,
                    task_structure=_task_structure)

    def add_task(self,task_name,task_type,task_structure):
        """
            Args:
                task_name(str):
                task_type(str):
                task_structure(str):
        """
        assert isinstance(task_name,str)
        assert isinstance(task_type,str)
        assert isinstance(task_structure,str)
        _base_directory = self.base_directory
        _task_name = task_name
        _task_type = task_type

        # <-------- determine task directory
        _task_directory = os.path.join(
                _base_directory,
                _task_name)

        # <-------- determine location of the structure
        _task_structure = os.path.join(
                self.structures['structure_directory'],
                self.structures['structures'][task_structure])

        # <-------- instantiate the task
        _module_name = TaskToClassMap[_task_type]['module']
        _class_name = TaskToClassMap[_task_type]['class']
        _module = importlib.import_module(_module_name)
        _class = getattr(_module,_class_name)
        self.obj_Task[_task_name] = _class(
                task_name = _task_name,
                task_directory = _task_directory,
                structure_filename = _task_structure,
                restart = False,
                fullauto = False)