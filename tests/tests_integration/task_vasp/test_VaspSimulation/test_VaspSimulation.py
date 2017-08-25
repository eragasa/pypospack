import os, shutil, subprocess
import ase.build
import pypospack.io.vasp as vasp
import pypospack.crystal as crystal
import pypospack.io.slurm as slurm
import pypospack.task.vasp as tsk_vasp

if __name__ == "__main__":
    # task information
    task_name = 'MgO_calc'
    task_directory = 'MgO_calc'
    # structure info
    structure_filename = "rsrc/MgO_NaCl_prim.vasp"
    print(os.path.exists(structure_filename))

    # changes to the standard run
    incar_config = {\
            'encut':800,
            'ismear':0,
            'sigma':0.05}
    
    # slurm info
    slurm_dict = {\
            'email':'eragasa@ufl.edu',
            'qos':'phillpot',
            'ntasks':16,
            'time':"1:00:00"}

    # initialize without path existing
    if os.path.exists(task_directory):
        shutil.rmtree(task_directory)

    task = tsk_vasp.VaspSimulation(\
            task_name='MgO_calc',
            task_directory='MgO_calc',
            restart=False)

    task.config(\
            poscar=structure_filename,
            incar=incar_config,
            xc='GGA')

    task.run(\
            job_type='slurm',
            exec_dict=slurm_dict)

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


