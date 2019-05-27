import os, shutil, subprocess
import ase.build
import pypospack.io.vasp as vasp
import pypospack.crystal as crystal
import pypospack.io.slurm as slurm
import pypospack.task.vasp as tsk_vasp

class TestVaspStructuralMinimization(object):

    def test_init(self):
        task_name = 'MgO_min'
        task_directory = 'MgO_min'
        structure_filename = 'rsrc/MgO_NaCl_prim.vasp'
        incar_config = {\
                'encut':800,
                'ismear':0,
                'sigma':0.05,
                'ispin':1}
        xc = 'GGA'

        task = tsk_vasp.VaspStructuralMinimization(\
                task_name=task_name,
                task_directory=task_directory,
                restart=False)

        assert os.path.exists(task_directory)
        assert task.task_name == task_name
        assert task.task_directory == os.path.join(os.getcwd(),task_directory)
        assert task.status == 'INIT'


    def test_config(self):
        task_name = 'MgO_min'
        task_directory = 'MgO_min'
        structure_filename = 'rsrc/MgO_NaCl_prim.vasp'
        incar_config = {\
                'encut':800,
                'ismear':0,
                'sigma':0.05,
                'ispin':1}
        xc = 'GGA'

        if os.path.exists(task_directory):
            shutil.rmtree(task_directory)

        task = tsk_vasp.VaspStructuralMinimization(\
                task_name=task_name,
                task_directory=task_directory,
                restart=False)
        task.config(poscar=structure_filename,incar=incar_config,xc=xc)

        assert task.incar.ibrion == 2
        assert task.incar.isif == 3

        assert os.path.exists(os.path.join(task.task_directory,'POSCAR'))
        assert os.path.exists(os.path.join(task.task_directory,'INCAR'))
        assert os.path.exists(os.path.join(task.task_directory,'POTCAR'))
        assert os.path.exists(os.path.join(task.task_directory,'KPOINTS'))

        poscar = vasp.Poscar();poscar.read(os.path.join(task.task_directory,'POSCAR'))
        incar = vasp.Incar();incar.read(os.path.join(task.task_directory,'INCAR'))
        potcar = vasp.Potcar();potcar.read(os.path.join(task.task_directory,'POTCAR'))
        kpoints = vasp.Kpoints();kpoints.read(os.path.join(task.task_directory,'KPOINTS'))


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
   # form of the exchange correlation functional 
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


