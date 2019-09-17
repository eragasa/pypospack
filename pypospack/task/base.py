import os,shutil
from collections import OrderedDict

class Task(object):
    """
    An abstract task for the implementation of pypospack.task objects.

    Args:
    Attributes:
        status_states(list of str): the list of status states
        task_name(str): the task name, will be used in
            job submission scripts
        task_directory(str): the directory the task,
            will be done
        status(str): The status states of a task is provided by
        conditions_INIT(OrderedDict)
        conditions_CONFIG(OrderedDict)
        conditions_READY(OrderedDict)
        conditions_RUNNING(OrderedDict)
        conditions_POST(OrderedDict)
        conditions_FINISHED(OrderedDict)
    Notes:


    """
    def __init__(self,
            task_name,
            task_directory,
            restart=False):
        # supported status states
        self.status_states = self.get_status_states()
        self.supported_status_states = self.get_status_states()
        # private member varaibles
        self._is_restart = None
        self._task_name = None
        self._task_directory = None
        self._status = None
        # process initialization arguments
        self.is_restart = restart
        self.task_name = task_name
        self.task_directory = task_directory
        self.root_directory = os.getcwd()
        self.results = None
        if restart is True:
            self.restart()
        else:
            self.create_task_directory(self.task_directory)
            self.update_status()

    def on_update_status(self):
        if self.status == 'INIT': 
            self.on_init()
        elif self.status == 'CONFIG': 
            self.on_config()
        elif self.status == 'READY': 
            self.on_ready()
        elif self.status == 'RUNNING':
            self.on_running()
        elif self.status == 'POST':
            self.on_post()
        elif self.status == "FINISHED":
            self.on_finished()
        elif self.status == "ERROR":
            self.on_error()
    
    def update_status(self):
        self.get_conditions_init()
        self.get_conditions_config()
        self.get_conditions_ready()
        self.get_conditions_running()
        self.get_conditions_post()
        self.get_conditions_finished()
        self.get_conditions_error()
        self.all_conditions_INIT \
                = all([v for k,v in self.conditions_INIT.items()])
        self.all_conditions_CONFIG \
                = all([v for k,v in self.conditions_CONFIG.items()])
        self.all_conditions_READY \
                = all([v for k,v in self.conditions_READY.items()])
        self.all_conditions_RUNNING \
                = all([v for k,v in self.conditions_RUNNING.items()])
        self.all_conditions_POST \
                = all([v for k,v in self.conditions_POST.items()])
        self.all_conditions_FINISHED \
                = all([v for k,v in self.conditions_FINISHED.items()])

        if self.all_conditions_INIT:
            self.status = 'INIT'
        if self.all_conditions_INIT and self.all_conditions_CONFIG:
            self.status = 'CONFIG'
        if self.all_conditions_INIT and self.all_conditions_CONFIG\
                and self.all_conditions_READY:
            self.status = 'READY'
        if self.all_conditions_INIT \
                and self.all_conditions_CONFIG\
                and self.all_conditions_READY \
                and self.all_conditions_RUNNING:
            self.status = 'RUNNING'
        if self.all_conditions_INIT \
                and self.all_conditions_CONFIG\
                and self.all_conditions_READY \
                and self.all_conditions_RUNNING\
                and self.all_conditions_POST:
            self.status = 'POST'
        if self.all_conditions_INIT \
                and self.all_conditions_CONFIG \
                and self.all_conditions_READY \
                and self.all_conditions_RUNNING \
                and self.all_conditions_POST \
                and self.all_conditions_FINISHED:
            self.status = 'FINISHED'
       
    def get_conditions_init(self):
        raise NotImplementedError

    def get_conditions_config(self):
        raise NotImplementedError

    def get_conditions_ready(self):
        raise NotImplementedError

    def get_conditions_running(self):
        raise NotImplementedError

    def get_conditions_post(self):
        raise NotImplementedError

    def get_conditions_finished(self):
        raise NotImplementedError

    def get_status_states(self):

        return ['INIT','CONFIG','READY','RUNNING','POST','FINISHED','ERROR']
   
    def on_init(self):
        raise NotImplementedError

    def on_config(self):
        raise NotImplementedError

    def on_ready(self):
        raise NotImplementedError

    def on_post(self):
        raise NotImplementedError

    def on_finished(self):
        raise NotImplementedError

    def on_error(self):
        raise NotImplementedError

    def restart(self):
        #<--- check if init has already occured
        if not os.path.exists(self.task_directory):
            self.create_task_directory(self.task_directory)

    def run(self):
        raise NotImplementedError()

    def create_task_directory(self,task_directory):
        #<--- check condition, we cannot run a simulation in the 
        #     same directory which we run this script.
        if os.path.abspath(self.root_directory) \
                == os.path.abspath(task_directory):
            err_msg = (
                "Cannot set the path of the task_directory to the root_directory.\n"
                "\ttask_directory:{}\n"
                "\troot_directory:{}\n").format(
                    task_directory,
                    root_directory)
            raise ValueError(err_msg)

        #<--- change to task_directory to an absolute path
        self.task_directory = os.path.abspath(task_directory) 
        
        #<--- we do this if the path already exists
        if os.path.isdir(self.task_directory):
            shutil.rmtree(self.task_directory,ignore_errors=True)
        
        #<--- hack to deal with race condition
        # https://stackoverflow.com/questions/12468022/python-fileexists-error-when-making-directory
        while True:
            try:
                os.mkdir(self.task_directory)
                break
            except FileExistsError as e:
                if os.path.isdir(self.task_directory):
                    shutil.rmtree(self.task_directory,ignore_errors=True)
                if e.errno != os.errno.EEXIST:
                    raise
                pass
        
        self.status = 'INIT'
    
    @property
    def task_name(self): 
        return self._task_name
    
    @task_name.setter
    def task_name(self, task_name):
        if isinstance(task_name,str):
            self._task_name = task_name
        else:
            self._task_name = str(task_name)

    @property
    def task_directory(self): 
        return self._task_directory

    @task_directory.setter
    def task_directory(self, task_directory):
        assert type(task_directory) is str
        self._task_directory = task_directory
    @property
    def is_restart(self):
        return self._is_restart

    @is_restart.setter
    def is_restart(self,is_restart):
        if not isinstance(is_restart,bool):
            msg_err = "is_restart must be a boolean value"
        self._is_restart = is_restart

    @property
    def status(self): 
        return self._status

    @status.setter
    def status(self,status):
        if status not in self.status_states:
            msg_err = (
                    "unsupported status state.\n"
                    "\tstatus:{}\n"
                    "\tsupported_status_states\n"
                    ).format(
                            status,
                            self.supported_status_states)
            raise ValueError(msg_err)
        self._status = status


