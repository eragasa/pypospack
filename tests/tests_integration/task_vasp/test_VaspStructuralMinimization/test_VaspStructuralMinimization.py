import os, shutil, subprocess
import ase.build
import pypospack.io.vasp as vasp
import pypospack.crystal as crystal
import pypospack.io.slurm as slurm
import pypospack.task.vasp as tsk_vasp

if __name__ == "__main__":
    # task information
    task_name = 'MgO_min'
    task_directory = 'MgO_min'
    # structure info
    structure_filename = "rsrc/MgO_NaCl_prim.vasp"

    # changes to the standard run
    incar_config = {\
            'encut':800,
            'ismear':0,
            'sigma':0.05,
            'ispin':1}
    

    xc = 'GGA'
    # slurm info
    slurm_dict = {\
            'email':'eragasa@ufl.edu',
            'qos':'phillpot',
            'ntasks':16,
            'time':"1:00:00"}

    task = tsk_vasp.VaspStructuralMinimization(\
            task_name='MgO_calc',
            task_directory='MgO_calc',
            restart=True)

    print('task.status:{}'.format(task.status))
    if task.status == 'INIT':
        task.config(poscar=structure_filename,incar=incar_config,xc=xc)

    if task.status == 'CONFIG':
        task.run(\
                job_type='slurm',
                exec_dict=slurm_dict)

    if task.status == 'RUN':
        pass

    if task.status == 'POST':
        task.postprocess()


    # commented out code
    if False:
        # initialize with path existing
        if os.path.exists(task_directory):
            shutil.rmtree(task_directory)
        os.mkdir(task_directory)
        task = task_vasp.VaspSimulation(\
                task_name='MgO_calc',
                task_directory='MgO_calc')


