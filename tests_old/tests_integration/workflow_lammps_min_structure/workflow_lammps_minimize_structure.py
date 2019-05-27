import os,sys
import importlib
import pypospack.workflow as workflow
import pypospack.task.lammps as lammps


if __name__ == "__main__":
    import inspect

    #for name, obj in inspect.getmembers(sys.modules[__name__]):
    #    print(name,object)

    wf_name = 'minimization'
    wf_directory = 'minimization'

    wf = workflow.WorkflowManager(wf_name=wf_name,wf_directory=wf_directory)
    print(__file__)
    wf.add_task('MgO_NaCl.E_min_all','lammps','LammpsSimulation')

    for k,v in wf.tasks.items():
        print(k,v['dir'])
