import unittest
import os, shutil, subprocess
import ase.build
import pypospack.io.vasp as vasp
import pypospack.crystal as crystal
import pypospack.io.slurm as slurm
import pypospack.task.vasp as tsk_vasp

class test_VaspSimulation_init(unittest.TestCase):

    def setUp(self):
        self.task_name = 'MgO_calc'
        self.task_directory = 'MgO_calc'
        self.structure_filename = 'rsrc/MgO_NaCl.prim.vasp'

        if os.path.exists(task_directory):
            shutil.rmtree(task_directory)

    def test_init_no_restart(self):
        task = tsk_vasp.VaspSimulation(\
                task_name=self.task_name,
                task_directory=self.task_directory,
                restart = False)

    def test_init_with_restart(self):
        task = tsk_vasp.VaspSimulation(\
                task_name=self.task_name,
                task_directory=self.task_directory,
                restart = True)

 
if __name__ == "__main__":
    # task information
    task_name = 'MgO_calc'
    task_directory = 'MgO_calc'
    # structure info
    structure_filename = "rsrc/MgO_NaCl_prim.vasp"

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
            restart=True)

    if task.status == 'INIT':
        task.config(\
                poscar=structure_filename,
                incar=incar_config,
                xc='GGA')

    if task.status == 'CONFIG':
        task.run(\
                job_type='slurm',
                exec_dict=slurm_dict)

    if task.status == 'RUN':
        # job is not complete yet
        print('job is not complete')

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


