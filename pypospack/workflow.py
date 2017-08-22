import os,sys,time,copy
import importlib
import inspect

class WorkflowManagerError(Exception):
    pass

class WorkflowManager(object):
    """

    Args:
       wf_name(str): workflow task name
       wf_directory(str): workflow directory name

    Attributes:
       wf
    """
    def __init__(self,wf_name,wf_directory=None):
        self.workflow_name = wf_name
        self.workflow_directory = None
        self.tasks= {}
        self.task_configs = {}
        self.task_requirements = {}
        self.task_status = {}
        self.param_dict = {}
        self.__set_workflow_directory(wf_directory=wf_directory)
        self.__log('wf_name:{}'.format(wf_name))

    def initialize_tasks(self,task_obj_config):
        """ initial tasks from a dictionary

        task_obj_config(:obj:`dict` of :obj:`dict`):  a dict of configuration
            dictionaries

            >>> task_obj_config = {}
            >>> task_obj_config['task1'] = {}
            >>> task_obj_config['task1']['module_name'] = 'pypospack.task.lammps'
            >>> task_obj_config['task1']['class_name'] = 'LammpsSimulation'
            >>> task_obj_config['task2'] = {}
            >>> ...
        """

        # reset the task attribute
        self.tasks = {}
        self.task_status = {}
        for task in self.tasks.keys():
            self.task_status[task] = None

        self.task_config = {}
        for task, task_config in task_obj_config.items():
            self.task_config[task] = {}
            for k,v in task_config.items():
                self.task_config[task][k] = v

        # iterate over the list of task and add them to the list
        # self.add_ask instantiates the class
        for task_name,task_info in task_obj_config.items():
            task_module_name = task_info['module_name']
            task_class_name = task_info['class_name']
            self.add_task(task_name,task_module_name,task_class_name)

    def run(self,wait_time=0.05):
        all_tasks_complete = False
        i = 1

        task_status_map = {
                'INIT':self.process_task_status_init,
                'CONFIG':self.process_task_status_config,
                'READY':self.process_task_status_ready,
                'RUN':self.process_task_status_run,
                'POST':self.process_task_status_post,
                'DONE':self.process_task_status_done }
        while not all_tasks_complete:
            # check to see if all the tasks are complete
            if all(obj['obj'].status == 'DONE' for k,obj in self.tasks.items()):
                all_tasks_complete = True
            # update status
            for task in self.tasks.keys():
                self.task_status[task] = self.tasks[task]['obj'].status
            # process status
            for task, status in self.task_status.items():
                task_status_map[status](task)

            time.sleep(wait_time)

    def process_task_status_init(self,task_name):
        print('INIT',task_name)
        config_list = self.tasks[task_name]['obj'].req_config()
        config_dict = {}
        for key in config_list:
            try:
                config_dict[key] = self.task_config[task_name][key]
            except KeyError as e:
                raise WorkflowManagerError(\
                    ('Task cannot be configured because config_dict is '
                    'missing necessary configuration information.\n'
                    '\ttask_name:{}\n'
                    '\ttask_config_field_name:{}\n'
                    ).format(task_name,str(e)))
        self.tasks[task_name]['obj'].send_config(config_dict)

    def process_task_status_config(self,task_name):
        print('CONFIG',task_name)

        # since we are at config, we want information to get to ready
        ready_list = self.tasks[task_name]['obj'].req_ready()
        ready_dict = {}
        for key in ready_list:
            if key == 'param_dict':
                ready_dict[key] = copy.deepcopy(self.param_dict)
            else:
                raise ValueError('{} need the information {} to get '
                        'to ready status.  WorkflowManager does not '
                        'have this information'.format(task_name,key))
        self.tasks[task_name]['obj'].send_ready(ready_dict)

    def process_task_status_ready(self,task_name):
        print('READY',task_name)
        run_list = self.tasks[task_name]['obj'].req_run()
        run_dict = {}
        for key in run_list:
            pass
        self.tasks[task_name]['obj'].send_run(run_dict)

    def process_task_status_run(self,task_name):
        print('RUN',task_name)
        post_list = self.tasks[task_name]['obj'].req_post()
        post_dict = {}
        for key in post_list:
            pass
        self.tasks[task_name]['obj'].send_post(post_dict)

    def process_task_status_post(self,task_name):
        print('POST',task_name)
        done_list = self.tasks[task_name]['obj'].req_done()
        done_dict = {}
        for key in done_list:
            pass
        self.tasks[task_name]['obj'].send_done(done_dict)

    def process_task_status_done(self,task_name):
        print('DONE',task_name)

    def process_task_status_unknown(self,task_name):
        raise WorkflowError('unknown Task.status')
    def add_task(self,task_name,task_module_name,task_class_name):
        """ add a task to the workflow

        Args:
            task_name(str): unique task_name id
            task_module_name(str): module name to load
            task_class_name(str): module class to load
        """

        task_directory = os.path.join(\
                self.workflow_directory,
                task_name)

        # tasks are instantiated dynamically
        # getattr(module,str_class_name)(task_name,task_directory)
        self.tasks[task_name] = {}
        self.tasks[task_name]['dir'] = task_directory
        try:
            task_module = importlib.import_module(task_module_name)
            task_class = getattr(task_module,task_class_name)
        except:
            raise

        self.tasks[task_name]['obj'] = task_class(
                task_name=task_name,
                task_directory=task_directory)

    def __log(self,msg):
        """ temporary logging method """
        print(msg)

    def __set_workflow_directory(self,wf_directory):
        if wf_directory is None:
            self.workflow_directory = os.path.join(\
                    os.getcwd(),
                    self.workflow_name)
        elif os.path.isabs(wf_directory):
            self.workflow_directory = os.path.join(\
                    wf_directory,
                    self.workflow_name)
        elif not os.path.isabs(wf_directory):
            self.workflow_directory = os.path.join(\
                    self.os.getcwd(),
                    wf_directory,
                    self.workflow_name)
        else:
            raise ValueError("do not know what to do with wf_directory:{}".format(wf_directory))

        # make workflow directory
        self.__log("setting the workfflow directory to {}".format(\
                self.workflow_directory))
        if not os.path.exists(self.workflow_directory):
            os.mkdir(self.workflow_directory)
        self.__log("wf_directory={}".format(self.workflow_directory))

class LammpsWorkflowManager(WorkflowManager):
    def __init__(self,wf_name,wf_directory):
        WorkflowManager.__init__(wf_name,wf_directory)
        self.status = 'INIT'
    def initialize_tasks(self,task_obj_config):
        super(LammpsWorkflowManager, self).initialize_tasks(task_obj_config)

