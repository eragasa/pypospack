import os,shutil,subprocess
import pypospack.io.vasp as vasp
import pypospack.io.slurm as slurm
from pypospack.task import Task

class VaspSimulationError():
    pass

class VaspSimulation(Task):
    def __init__(self,task_name,task_directory,restart=True):
        Task.__init__(self,task_name,task_directory,restart)

        self.poscar = vasp.Poscar()
        self.incar = vasp.Incar()
        self.potcar = vasp.Potcar()
        self.kpoints = vasp.Kpoints()

        # additional items to add
        additional_config_dict = {\
                'incar':self.set_incar,
                'poscar':self.set_poscar,
                'xc':self.set_xc,
                'encut':self.set_encut}
        additional_ready_dict = {}
        additional_run_dict = {}
        additional_post_dict = {}

        # update configuration dictionaries
        self.config_dict.update(additional_config_dict)
        self.ready_dict.update(additional_config_dict)
        self.run_dict.update(additional_ready_dict)
        self.post_dict.update(additional_post_dict)

        self.status = 'INIT'

    def restart(self):
        """ overwrite original method """
        poscar_exists = os.path.exists(os.path.join(\
                self.task_directory,'POSCAR'))
        incar_exists = os.path.exists(os.path.join(\
                self.task_directory,'INCAR'))
        kpoints_exists = os.path.exists(os.path.join(\
                self.task_directory,'KPOINTS'))
        potcar_exists = os.path.exists(os.path.join(\
                self.task_directory,'POTCAR'))
        if poscar_exists and incar_exists and kpoints_exists and potcar_exists:
            self.config = 'CONFIG'
        else:
            shutil.rmtree(self.task_directory)
            os.mkdir(self.task_directory)
            self.config = 'INIT'

        if os.path.exists(os.path.join(\
                self.task_directory,'jobSubmitted')):
            self.config = 'RUN'

        if os.path.exists(os.path.join(\
                self.task_directory,'jobComplete')):
            self.config = 'POST'


    def config(self,poscar=None,incar=None,kpoints=None,xc='GGA'):
        # read the poscar file, then write it
        if poscar is not None:
            if isinstance(poscar,str):
                self.poscar_filename = poscar
                self.poscar.read(poscar)
            elif isinstance(poscar,pypospack.crystal.SimulationCell):
                self.poscar = Poscar(poscar)
            else:
                msg = (\
                        "argument 'poscar' must be a filename or an instance "
                        "of pypospack.crystal.SimulationCell, passed in {}"
                      ).format(poscar)
                raise KeyError
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
        if isinstance(incar,str):
            self.incar.read(incar)
        elif isinstance(incar,dict):
            for k,v in incar.items():
                setattr(self,'incar.{}'.format(k),v)
        elif isinstance(incar,pypospack.io.vasp.Incar):
            self.incar = vasp.Incar(incar)
        elif incar is None:
            pass
        else:
            raise VaspSimulationError

        if self.incar.encut == 'Auto':
            self.incar.encut = self.get_encut_from_potcar(xc)

        try:
            self.incar.write(\
                os.path.join(\
                    self.task_directory,'INCAR'))
        except:
            raise

        # process KPOINTS file
        if kpoints is None:
            pass
        elif isinstance(kpoints,str):
            self.kpoints.read(kpoints)
        elif isinstance(kpoints,dict):
            setattr(self,'kpoints.{}'.format(k),v)
        elif isinstance(kpoints,pypospack.io.vasp.Kpoints):
            self.kpoints = vasp.Kpoints(kpoints)
        else:
            raise VaspSimulationError

        self.kpoints.write(\
                os.path.join(\
                    self.task_directory,'KPOINTS'))

        self.status = 'CONFIG'

    def run(self,job_type='slurm',exec_dict=[]):
        self.job_type = job_type
        if self.job_type == 'slurm':
            fname = os.path.join(\
                    self.task_directory,'jobSubmitted')
            with open(fname,'w') as f:
                f.write('job started')

            slurm.write_vasp_batch_script(
                    filename=os.path.join(\
                            self.task_directory,
                            'runjob_hpg.slurm'),
                    job_name=self.task_name,
                    email=exec_dict['email'],
                    qos=exec_dict['qos'],
                    ntasks=exec_dict['ntasks'],
                    time=exec_dict['time'])

            cmd = 'sbatch runjob_hpg.slurm'

            old_wd = os.getcwd()
            os.chdir(self.task_directory)
            subprocess.call(cmd,shell=True)
            os.chdir(old_wd)
        else:
            msg = ('do not how to run the simulation using job_type'
                   '={}').format(job_type)
            raise VaspSimulationError

        self.status = 'RUN'

    def postprocess(self):
        pass

    def slurm_is_done(self):
        if os.path.exists(os.path.join(self.task_dir,'jobComplete')):
            self.status = 'POST'

    def set_incar(self,incar_info):
        cond1 = isinstance(incar_info,str)
        cond2 = os.path.isabs(incar_info)
        if cond1 and cond2:
            self.incar.read(incar)
        else:
            raise ValueError

    def set_poscar(self,poscar_info):
        cond1 = isinstance(incar_info,str)
        cond2 = os.path.isabs(incar_info)

        if cond1 and cond2:
            self.poscar.read(poscar)
        elif isinstance(poscar,pypospack.crystal.SimulationCell):
            self.poscar = vasp.Poscar(poscar)
        else:
            raise ValueError

    def set_xc(self,xc):
        self.xc = xc

    def set_encut(self,encut):
        self.encut = encut
    def req_config(self):
        pass

    def snd_config(self,config_dict):
        """ sends configuration information and configures
        """

        pass

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



