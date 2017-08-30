import os, shutil, subprocess
import ase.build
import pypospack.io.vasp as vasp
import pypospack.crystal as crystal
import pypospack.io.slurm as slurm
import pypospack.task.vasp as tsk_vasp

class TestVaspSimulation(object):

    def init_wo_restart(self):
        print('starting vasp simulation without restart')
        self.task_name = 'no_restart'
        self.task_directory = 'no_restart'
        # remove existing directory if it exists
        if os.path.exists(self.task_directory):
            shutil.rmtree(self.task_directory)

        self.structure_filename = "rsrc/MgO_NaCl_prim.vasp"
        self.task = tsk_vasp.VaspSimulation(
            task_name=self.task_name,
            task_directory=self.task_directory,
            restart=False)


    def test_init_wo_restart(self):
        self.init_wo_restart()
        assert self.task.task_name == self.task_name
        assert self.task.task_directory == os.path.join(os.getcwd(),self.task_directory)
        assert os.path.exists(self.task_directory)
        assert isinstance(self.task.poscar,vasp.Poscar)
        assert isinstance(self.task.incar,vasp.Incar)
        assert isinstance(self.task.potcar,vasp.Potcar)
        assert isinstance(self.task.kpoints,vasp.Kpoints)

    def test_config_no_param(self):
        self.init_wo_restart()
        self.task.config(poscar=self.structure_filename)

    def test_config_xc_gga(self):
        self.init_wo_restart()
        self.task.config(poscar=self.structure_filename,xc='GGA')

    def test_config_xc_lda(self):
        self.init_wo_restart()
        self.task.config(poscar=self.structure_filename,xc='LDA')

    def test_config_change_encut(self):
        self.init_wo_restart()
        incar_config = {'encut':123.}
        self.task.config(\
                poscar=self.structure_filename,
                incar=incar_config)
        assert self.task.incar.encut == 123.
while False:

    def test_config_change_incar(self):
        task_name = 'MgO_calc'
        task_directory = 'MgO_calc'
        if os.path.exists(task_directory):
            shutil.rmtree(task_directory)

        structure_filename = "rsrc/MgO_NaCl_prim.vasp"
        task = tsk_vasp.VaspSimulation(
            task_name='MgO_calc',
            task_directory='MgO_calc',
            restart=False)

        xc = 'GGA'
        task.config(poscar=structure_filename,incar=incar_config,xc=xc)
 
if __name__ == "__main__":

    # changes to the standard run
    incar_config = {\
            'encut':800,
            'ismear':0,
            'sigma':0.05}
    

    xc = 'GGA'
    # slurm info
    slurm_dict = {\
            'email':'eragasa@ufl.edu',
            'qos':'phillpot',
            'ntasks':16,
            'time':"1:00:00"}

    task = tsk_vasp.VaspSimulation(\
            task_name='MgO_calc',
            task_directory='MgO_calc',
            restart=True)

    print('task.status:{}'.format(task.status))

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


