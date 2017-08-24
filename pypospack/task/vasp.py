import pypospack.io.vasp as io_vasp
from pypospack.task import Task

class VaspSimulationError():
    pass

class VaspSimulation(Task):
    def __init__(self,task_name,task_directory,restart=True):
        Task.__init__(self,task_name,task_directory,restart)
        self.poscar = io_vasp.Poscar()
        self.incar = io_vasp.Incar()
        self.potcar = io_vasp.Potcar()
        self.kpoints = vasp.Kpoints()
        self.status = 'INIT'

        # adding additional items
        config_dict = {\
                'incar',self.set_incar
                'poscar',self.set_poscar
                'xc':self.set_xc,
                'encut':self.set_encut}
        ready_dict = {}
        run_dict = {}
        post_dict = {}

        self.config_dict = dict(self.config_dict.items()+config_dict.items()}
    def req_config(self):
        pass

    def snd_config(self,config_dict):
        """ sends configuration information and configures
        """

        pass

    def config(self,poscar,incar=None,xc='GGA',encut=None):
        # read the poscar file, then write it
        self.poscar_filename = pioscar
        self.poscar.read(poscar)
        self.poscar.write(\
                os.path.join(\
                    self.task_directory,'POSCAR'))

        # write the potcar file, the read in information
        self.potcar.symbols = self.poscar.symbols
        self.potcar.xc = xc
        self.potcar.write(\
                os.path.join(\
                    self.task_directory,"POTCAR"))
        self.potcar.read(\
                os.path.join(\
                    self.task_directory,"POTCAR"))

        # process incar file
        if incar is not None:
            self.incar.read(incar)
        if encut is not None:
            if encut == 'Auto':
                self.incar.encut = self.get_encut_from_potcar(xc)
            else:
                self.incar.encut = encut 
            self.incar.encut = encut
        try:
            self.incar.write(\
                os.path.join(\
                    self.task_directory,'INCAR'))
        except:
            raise

        # process KPOINTS file
        self.kpoints.write(\
                os.path.join(\
                    self.kpoints,'KPOINTS'))

    def run(self,job_type='slurm',exec_dict=[]):
        if job_type == 'slurm':
            fname = os.path.join(\
                    self.task_directory,'jobSubmitted')
            with open(fname,'w') as f:
                f.write('job started')

            slurm.write_vasp_batch_script(
                    filename=os.path.join(\
                            self.task_directory,
                            'runjo_hpg.slurm'),
                    job_name=self.task_name,
                    email=exec_dict[email]
                    qos=exec_dict[qos]
                    ntasks=exec_dict[qos]
                    time=exec_dict[time])

            cmd = 'sbatch runjob_hpg.slurm'
            subprocess.call(cmd,shell=True)
        else:
            msg = 'do not how to run the simulation using job_type={}'.format(job_type)
            raise VaspSimulationError

    def set_xc(self,xc):
        self.xc = xc

    def set_encut(self,encut):
        self.encut = encut

    def set_incar(self,incar):
        if isinstance(incar,str):
            self.incar.read(incar)
        elif isinstance(incar,dict):
            for k in incar:
                setattr(self.incar,k,incar[k])
        else:
            msg = 'do not know how to process incar params'
            raise VaspSimulationError(msg)
    def get_encut_tfrom_potcar(self):
        return max(self.potcar.encut_max)



