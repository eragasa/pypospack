import os, shutil

class Workflow():

    def get_workflow_name(structure_name):
        msg = "subclasses of workflow need to implement this method"
        raise NotImplementedError(msg)

    def create_tasks(self):
        msg = "subclasses of workflow need to implement this method"
        raise NotImplementedError(msg)

    def create_workflow_path(self,path):
        if os.path.isdir(path):
            assert path is not os.getcwd()
            shutil.rmtree(path)
        os.mkdir(path)
