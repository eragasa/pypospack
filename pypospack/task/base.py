class Task(object):
    """
    Attributes:
        task_name(str): the task name, will be used in
            job submission scripts
        task_directory(str): the directory the task,
            will be done
    Notes: 
    If the __init__ is overridden, the attributes
    task_name, task_directory, and status must be 
    set.

    """
    def __init__(self,task_name,task_directory,restart=False):
        self.task_name = task_name
        self.task_directory = None

        self.config_dict = {}
        self.ready_dict = {}
        self.run_dict = {}
        self.post_dict = {}
        if os.path.isabs(task_directory):
            # absolute path, set directly.
            self.task_directory = task_directory
        else:
            # if a relative path, make it an
            # absolute path
            self.task_directory = os.path.join(\
                    os.getcwd(),
                    task_directory)

        # check to see if path exists
        if os.path.exists(self.task_directory):
            if restart:
                self.restart()
            else:
                # if no restart, start
                shutil.rmtree(self.task_directory)
        os.mkdir(self.task_directory)
    def restart(self):
        raise NotImplementedError


        self.task_directory = os.path.join(os.getcwdtask_directory
        self.status = 'INIT'

class TestTask(object):
