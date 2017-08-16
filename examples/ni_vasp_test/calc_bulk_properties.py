import pypospack.io.vasp as vasp

class PypospackWorkflow(object):
    def __init__(self):
        self.tasklist = {}

    def add_task(self,taskname,obj_task):
        self.tasklist[taskname] = obj_task

class PypospackTask(object):
    def __init__(self):
        self.status = 'NOTSTARTED'
        self.taskname = None
    def add_child(self,obj):
        assert isinstance(obj,PypospackTask)

class CalculateBulkProperties(PypospackTask):
    def __init__(self,taskname=None,structure=None):
        PypospackTask.__init__(self)
        self.taskname = taskname
        self.p_structure = structure
    
    def config(self):
        self.status = 'CONFIG'

    def run(self):

        self.status = 'WAIT'        

    def postprocess(self):
        self.status = 'COMPLETE'

class VaspMinimizeStructure(PypospackTask):
    def __init__(self):
        Pypospack.__init__(self,param_structure=None)
        self.vaspsimulation = vasp.VaspSimulation()
        self.simulation_directory = self.task_dir
        if param_structure is not None:
        self.vaspsimulation.poscar = vasp.Poscar(\
                ase.build.bulk(param_structure[0],
                               param_structure[1],
                               param_structure[2]
                               param_structure[3]
         

if __name__ == '__main__':

 
    structures = {}
    structures['Ni_fcc_cubic'] = [['Ni'],'fcc', 3.508,'cubic']
    structures['Ni_bcc_cubic'] = [['Ni'],'bcc', 3.508,'cubic']
    structures['Ni_hcp_cubic'] = [['Ni'],'hcp', 3.508,'cubic']
    structures['Ni_dia_cubic'] = [['Ni'],'diamond', 3.508,'cubic']
    structures['Ni_sc_cubic']  = [['Ni'],'sc', 3.508,'cubic']
    
    wf = PypospackWorkflow()
    for k,v in structures.items():
        wf.add_task(k,CalculateBulkProperties(k,structures[k]))
        wf.tasklist[k].config(structures[k]) 

