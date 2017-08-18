import os,sys
import importlib
import pypospack.workflow as workflow
import pypospack.task as task 

#class Task(object):
#    def __init__(self,task_name,task_directory):
#        self.task_name = task_name
#        self.task_directory = task_directory
#
#    def request_requirements(self):
#        return NotImplementedError

class DerivedTask(pypospack.task.Task):
    def __init__(self,task_name,task_directory):
        Task.__init__(self,task_name,task_directory)

if __name__ == "__main__":
    import inspect

    #for name, obj in inspect.getmembers(sys.modules[__name__]):
    #    print(name,object)

    wf_name = 'wf_test'
    wf_directory = 'wf_test'
    wf = workflow.WorkflowManager(wf_name=wf_name,wf_directory=None)
    print(__file__)
    wf.add_task('derived_task','__main__','DerivedTask')
    wf.add_task('derived_task2','__main__','DerivedTask')

    for k,v in wf.tasks.items():
        print(k,v['dir'])
