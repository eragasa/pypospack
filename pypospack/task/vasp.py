import os,shutil,subprocess
import pypospack.io.vasp as vasp
import pypospack.io.slurm as slurm
from pypospack.task import Task

class VaspSimulationError():
    pass

class VaspSimulation(Task):
    """

    The default configuration will do calculate the Energy of the structure
    using the GGA functional.  The energy cutoff for the plane waves will be
    set at the highest ENMAX contained in the potcar files.  The Kpoints will
    be set at 6x6x6 kpoint mesh in the Brillioun zone.

    Args:
        task_name(str):
        task_directory(str): where to create a directory
        restart(bool): if this flag is set to True, then it will
            attempt to restart the simulation from information 
            contained in task_directory.  Default is True.

    Attributes:
        poscar(pypospack.io.vasp.Poscar): A class to manage poscar files
        incar(pypospack.io.vasp.Incar): A class to manage incar files
        potcar(pypospack.io.Potcar): A class to manage potcar files
        kpoints(pypospak.io.Kpoints): A class to manage kpoints

    """
    def __init__(self,task_name,task_directory,restart=True):
        Task.__init__(self,task_name,task_directory,restart)

        if self.status == 'INIT':
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

        print('restarting simulation {}'.format(self.task_name))
        if poscar_exists and incar_exists and kpoints_exists and potcar_exists:
            print('\tsimulation already configured...')
            self.status = 'CONFIG'
        else:
            print('\tsimulation not configured, initializing')
            shutil.rmtree(self.task_directory)
            os.mkdir(self.task_directory)
            self.status = 'INIT'

        if os.path.exists(os.path.join(\
                self.task_directory,'jobSubmitted')):
            print('\tjob has been submitted')
            self.status = 'RUN'
        else:
            print('\tjob has not been submitted')

        if os.path.exists(os.path.join(\
                self.task_directory,'jobCompleted')):
            print('\tjob is done')
            self.status = 'POST'
        else:
            print('job has not finished')

    def config(self,poscar=None,incar=None,kpoints=None,xc='GGA'):
        # read the poscar file, then write it
        self.read_poscar(poscar=poscar)
        self.write_poscar()

        # write the potcar file, the read in information
        self.write_potcar(xc=xc)
        self.read_potcar()

        # process incar file
        self.read_incar(incar=incar)
        self.config_incar()
        self.write_incar()
        # process KPOINTS file
        self.read_kpoints(kpoints)
        self.write_kpoints()

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

    def read_poscar(self,poscar=None):
        if poscar is not None:
            if isinstance(poscar,str):
                self.poscar_filename = poscar
                self.poscar.read(poscar)
            elif isinstance(poscar,pypospack.crystal.SimulationCell):
                self.poscar = vasp.Poscar(poscar)
            else:
                msg = (\
                        "argument 'poscar' must be a filename or an instance "
                        "of pypospack.crystal.SimulationCell, passed in {}"
                      ).format(poscar)
                raise KeyError

    def write_poscar(self):
        self.poscar.write(\
                os.path.join(\
                    self.task_directory,'POSCAR'))

    def read_potcar(self):
        self.potcar.read(\
                os.path.join(\
                    self.task_directory,"POTCAR"))

    def write_potcar(self,xc='GGA'):
        self.potcar.symbols = self.poscar.symbols
        self.potcar.xc = xc
        self.potcar.write(\
                os.path.join(\
                    self.task_directory,"POTCAR"))

    def read_incar(self,incar='None'):
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

    def config_incar(self):
        potcar_encut_min = max(self.potcar.encut_min)
        potcar_encut_max = 2*max(self.potcar.encut_max)
        if self.incar.encut == 'Auto':
            self.incar.encut = max(self.potcar.encut_max)
        if self.incar.encut >= potcar_encut_min:
            self.incar.encut = potcar_encut_min
        if self.incar.encut <= potcar_encut_max:
            self.incar.encut = potcar_encut_max

    def write_incar(self):
        try:
            self.incar.write(\
                os.path.join(\
                    self.task_directory,'INCAR'))
        except:
            raise

    def read_kpoints(self,kpoints=None):
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

    def config_kpoints(self):

        # check for tetrahedron method and switch to
        # gamma
        if self.incar.ismear in [-4,-5]:
            self.kpoints.mesh_type = 'Gamma'

        # check to see if we have a hexagonal cell
        
    def write_kpoints(self):
        self.kpoints.write(\
                os.path.join(\
                    self.task_directory,'KPOINTS'))
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


class VaspStructuralMinimization(VaspSimulation):
    def __init__(self,task_name,task_directory,restart=True):
        VaspSimulation.__init__(self,task_name,task_directory,restart)

    def config_incar(self):
        VaspSimulation.config_incar(self)

        if self.incar.ibrion not in [1,2]:
            self.ibrion = 2 # set to conjugate gradient method

        self.isif = 3 # relax everything

