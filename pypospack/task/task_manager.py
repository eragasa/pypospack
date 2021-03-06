import os,signal,copy,importlib,time
from collections import OrderedDict

from pypospack.exceptions import PypospackTaskManagerError

TaskToClassMap = OrderedDict()
TaskToClassMap['lmps_min_all'] = OrderedDict()
TaskToClassMap['lmps_min_all']['module'] = 'pypospack.task.lammps'
TaskToClassMap['lmps_min_all']['class'] = 'LammpsStructuralMinimization' 
TaskToClassMap['lmps_min_pos'] = OrderedDict()
TaskToClassMap['lmps_min_pos']['module'] = 'pypospack.task.lammps'
TaskToClassMap['lmps_min_pos']['class'] = 'LammpsPositionMinimization'
TaskToClassMap['lmps_min_none'] = OrderedDict()
TaskToClassMap['lmps_min_none']['module'] = 'pypospack.task.lammps'
TaskToClassMap['lmps_min_none']['class'] = 'LammpsStaticCalculations'
TaskToClassMap['lmps_min_sf'] = OrderedDict()
TaskToClassMap['lmps_min_sf']['module'] = 'pypospack.task.lammps'
TaskToClassMap['lmps_min_sf']['class'] = 'LammpsStackingFaultMinimization'
TaskToClassMap['lmps_elastic'] = OrderedDict()
TaskToClassMap['lmps_elastic']['module'] = 'pypospack.task.lammps'
TaskToClassMap['lmps_elastic']['class'] = 'LammpsElasticCalculation'
TaskToClassMap['lmps_npt'] = OrderedDict()
TaskToClassMap['lmps_npt']['module'] = 'pypospack.task.lammps'
TaskToClassMap['lmps_npt']['class']= 'LammpsNptSimulation'
TaskToClassMap['lmps_neb'] = OrderedDict()
TaskToClassMap['lmps_neb']['module'] = 'pypospack.task.lammps'
TaskToClassMap['lmps_neb']['class'] = 'LammpsNebCalculation'
TaskToClassMap['gulp_gamma_phonons'] = OrderedDict()
TaskToClassMap['gulp_gamma_phonons']['module'] = 'pypospack.task.gulp'
TaskToClassMap['gulp_gamma_phonons']['class'] = 'GulpGammaPointPhonons'


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

    def __all_simulations_finished(self,obj_Tasks):
        """

        Returns:
            True if all simulations are completed
            False if all simulations are not completed
        """
        _statuses = [v.status in ['FINISHED','ERROR'] for v in obj_Tasks.values()]
        return all(_statuses) 

    def make_eam_setfl_file(self,potential,filename):
        pass


    def evaluate_tasks(self,parameters,potential,max_time_per_simulation=100):
        """
        parameters (OrderedDict) key-value pair of parameter names and parameter values
        potential (OrderedDict) key-value pair of potential descriptors
        max_time_per_simulation (int) 
        """
        
        _sleep_time = 0.1
        _max_time_per_simulation = max_time_per_simulation

        # initialize results dictions
        self.results = OrderedDict()

        # each task requires potential information and parameter information provided in a dictionary
        _configuration = OrderedDict()
        _configuration['potential'] = potential
        _configuration['parameters'] = parameters
  
        _start_time = time.time()
        while not self.__all_simulations_finished(self.obj_Task):
            
            # if the maximum time has been exceeded for this parameter set, we are going to kill
            # off all the subprocesses which maybe running simulations in each of the tasks.
            _time_elapsed = time.time() - _start_time
            if _time_elapsed > _max_time_per_simulation:
                for k_task,o_task in self.obj_Task.items():
                    # kill off process
                    # https://www.programcreek.com/python/example/11892/os.getpgid
                    # https://stackoverflow.com/questions/4789837/how-to-terminate-a-python-subprocess-launched-with-shell-true/4791612#4791612
                    # https://www.codeday.top/2017/06/28/25301.html
                    try:
                        o_task.process.kill()
                        #pid = o_task.process.pid
                        #pgid = os.getpgid(pid)
                        #if pgid == pid:
                        #    os.killpg(pgid,signal.SIGTERM)
                        #else:
                        #    os.kill(pgid,signal.SIGTERM)
                    except: 
                        pass
                raise PypospackTaskManagerError('simulation time exceeded',parameters=parameters)
            
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

                    _configuration = OrderedDict()
                    _configuration['potential'] = potential
                    _configuration['parameters'] = parameters
                    if 'bulk_structure' in self.tasks[k_task]:
                        _structure_name = self.tasks[k_task]['bulk_structure']
                        _structure_filename = os.path.join(
                                self.structures['structure_directory'],
                                self.structures['structures'][_structure_name])
                        _configuration['bulk_structure'] = _structure_name
                        _configuration['bulk_structure_filename'] = _structure_filename
                    
                    o_task.on_init(configuration=_configuration)

                elif o_task.status == 'CONFIG':
                    try:
                        o_task.on_config(
                            configuration=_configuration,
                            results=self.results)
                    except TypeError as e:
                        o_task.on_config(configuration=_configuration)
                elif o_task.status == 'READY':
                    try:
                        o_task.on_ready(results=self.results)
                    except TypeError as e:
                        print("Error with {}:{}".format(k_task,type(o_task)))
                        raise
                elif o_task.status == 'RUNNING':
                    o_task.on_running()
                elif o_task.status == 'POST':
                    o_task.on_post()
                    _results = o_task.results
                    try:
                        for k,v in o_task.results.items():
                            self.results[k] = v
                    except AttributeError as e:
                        print('k_task:{}'.format(k_task))
                        print('o_task:{}'.format(o_task))
                        raise

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
            try:
                _task_options = v_task['task_options']
            except KeyError as e:
                _task_options = None

            self.add_task(
                    task_name=_task_name,
                    task_type=_task_type,
                    task_structure=_task_structure,
                    task_options = _task_options)

    def add_task(self,task_name,task_type,task_structure,task_options=None):
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
        _task_options = copy.deepcopy(task_options)
        # <-------- determine task directory
        _task_directory = os.path.join(
                _base_directory,
                _task_name)

        # <-------- determine location of the structure
        _structure_dir = self.structures['structure_directory']
        
        try:
            _structure_fn = self.structures['structures'][task_structure]
        except KeyError as e:
            print('the structure name \'{}\' is not defined'.format(task_structure))
            raise

        _task_structure = os.path.join(_structure_dir,_structure_fn)

        # <-------- instantiate the task
        _module_name = TaskToClassMap[_task_type]['module']
        _class_name = TaskToClassMap[_task_type]['class']
        _module = importlib.import_module(_module_name)
        _class = getattr(_module,_class_name)
        
        if _task_options is None:
            self.obj_Task[_task_name] = _class(
                    task_name = _task_name,
                    task_directory = _task_directory,
                    structure_filename = _task_structure,
                    restart = False,
                    fullauto = False)
        else:
            kwargs = _task_options
            self.obj_Task[_task_name] = _class(
                    task_name = _task_name,
                    task_directory = _task_directory,
                    structure_filename = _task_structure,
                    restart = False,
                    fullauto = False,
                    **kwargs)
