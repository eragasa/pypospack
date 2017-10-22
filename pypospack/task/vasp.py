# -*- coding: utf-8 -*-
"""Input and output functions and classes for VASP """
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2016,2017"
__license__ = "Simplified BSD License"
__version__ = "1.0"

import os,shutil,subprocess,yaml,copy,pathlib
import numpy as np
import pypospack.crystal as crystal
import pypospack.io.vasp as vasp
import pypospack.io.slurm as slurm
from pypospack.task import Task

class VaspSimulationError(Exception):
    """Error class for dealing with VASP simulation issues"""
    def __init_(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)

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

        self.poscar = vasp.Poscar()
        self.incar = vasp.Incar()
        self.potcar = vasp.Potcar()
        self.kpoints = vasp.Kpoints()
        if self.status == 'INIT':

            # additional items to add
            additional_config_dict = {\
                'incar':self.set_incar,
                'poscar':self.set_poscar,
                'kpoints':self.set_kpoints,
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

        if os.path.exists(os.path.join(self.task_directory,'jobPost')):
            self.status = 'DONE'
        elif os.path.exists(os.path.join(self.task_directory,'jobCompleted')):
            self.status = 'POST'
        elif os.path.exists(os.path.join(self.task_directory,'jobSubmitted')):
            self.status = 'RUN'
        elif poscar_exists and incar_exists and kpoints_exists and potcar_exists:
            self.status = 'CONFIG'
        else:
            shutil.rmtree(self.task_directory)
            os.mkdir(self.task_directory)
            self.status = 'INIT'

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
        """Some postprocess steps

        - removes the POTCAR file, for copyright reasons
        - reads in the OUTCAR file
        - reads in the OSZICAR file
        - reads in the CONTCAR file

        """
        try:
            os.remove(os.path.join(self.task_directory,'POTCAR'))
        except FileNotFoundError:
            pass
        self.outcar = vasp.Outcar()
        self.outcar.read(os.path.join(self.task_directory,'OUTCAR'))
        self.contcar = vasp.Poscar()
        self.contcar.read(os.path.join(self.task_directory,'CONTCAR'))
        # self.oszicar = vasp.Oszicar()
        # self.oszicar.read(os.path.join(self.task_directory,'OSZICAR'))

    def read_poscar(self,poscar=None):
        """ read in a poscar file or object

        This method reads in new structural information into the class.  It does 
        not write the POSCAR file which may already exist in the the simulation
        directory.

        Args:
            poscar (:obj:`str` or :obj:`pypospack.crystal.SimulationCell`): 
                structure information file either provided as the location of a 
                POSCAR file, or from an instance of 
                :obj:`pypospack.crystal.SimulationCell`

        """
        if poscar is None:
            poscar = os.path.join(self.task_directory,'POSCAR')
        
        if isinstance(poscar,crystal.SimulationCell):
            self.poscar = vasp.Poscar(poscar)
        elif isinstance(poscar,str):
            self.poscar = vasp.Poscar()
            self.poscar.read(poscar)
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
        if incar is None:
            if os.path.exists(os.path.join(self.task_directory,'INCAR')):
                self.incar = vasp.Incar()
                self.incar.read(os.path.join(self.task_directory,'INCAR'))
            else:
                self.incar = vasp.Incar()
        elif isinstance(incar,str):
            self.incar = vasp.Incar()
            self.incar.read(incar)
        elif isinstance(incar,dict):
            # check to see if the incar already exists, and read it in
            if os.path.exists(os.path.join(self.task_directory,'INCAR')):
                self.incar = vasp.Incar()
                self.incar.read(os.path.join(self.task_directory,'INCAR'))
            # otherwise, we're just going to use the default configuarion
            else:
                self.incar = vasp.Incar()
            # now update based upon what is dictionary
            for k,v in incar.items():
                setattr(self.incar,k,v)
        else:
            raise VaspSimulationError

    def config_incar(self):
        pass
        #potcar_encut_min = max(self.potcar.encut_min)
        #potcar_encut_max = 2*max(self.potcar.encut_max)
        #if self.incar.encut == 'Auto':
        #    self.incar.encut = max(self.potcar.encut_max)
        #if self.incar.encut >= potcar_encut_min:
        #    self.incar.encut = potcar_encut_min
        #if self.incar.encut <= potcar_encut_max:
        #    self.incar.encut = potcar_encut_max

    def write_incar(self):
        try:
            self.incar.write(\
                os.path.join(\
                    self.task_directory,'INCAR'))
        except:
            raise

    def read_kpoints(self,kpoints=None):
        if kpoints is None:
            if os.path.exists(os.path.join(self.task_directory,'KPOINTS')):
                self.kpoints = vasp.Kpoints()
                self.kpoints.read(os.path.join(self.task_directory,'KPOINTS'))
            else:
                self.kpoints = vasp.Kpoints()
        elif isinstance(kpoints,str):
            self.kpoints = vasp.Kpoints()
            self.kpoints.read(kpoints)
        elif isinstance(kpoints,dict):
            if os.path.exists(os.path.join(self.task_directory,'KPOINTS')):
                self.kpoints = vasp.Kpoints()
                self.kpoints.read(os.path.join(self.task_directory,'KPOINTS'))
            else:
                self.kpoints = vasp.Kpoints()
            for k,v in kpoints.items():
                setattr(self.kpoints,k,v)
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

    def set_kpoints(self,kpoints_info):
        raise NotImplementedError

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
    def get_encut_from_potcar(self):
        return max(self.potcar.encut_max)

class VaspStructuralMinimization(VaspSimulation):
    """ Task class for VASP structural minimization

    The default configuration will do a full structural minimization using
    the GGA functional.  The energy cutoff for the plane waves will be set
    at the highest ENMAX contained in the POTCAR files.  The Kpoints will be
    set at the 6x6x6 kpoint mesh in the Brillioun zone.  Forces are 
    minimized to 1e-3 ev/Angs.

    Args:
        task_name(str): the task name for this VASP simulation
        task_directory(str): where to create a directory
        restart(bool): if this flag is set to True, then it will attempt
           restart the simulation from information contained in task_directory.
           Default is True.

    Attributes:
        poscar(pypospack.io.vasp.Poscar): A class to manage poscar files
        incar(pypospack.io.vasp.Incar): A class to manage incar files
        potcar(pypospack.io.Potcar): A class to manage potcar files
        kpoints(pypospak.io.Kpoints): A class to manage kpoints

    """
    def __init__(self,task_name,task_directory,restart=True):
        VaspSimulation.__init__(self,task_name,task_directory,restart)

    def config_incar(self):
        VaspSimulation.config_incar(self)

        if self.incar.ibrion not in [1,2]:
            self.incar.ibrion = 2 # set to conjugate gradient method
        self.incar.isif = 3 # relax everything
        self.incar.write(os.path.join(self.task_directory,'INCAR'))

class VaspPositionMinimization(VaspSimulation):
    def __init__(self,task_name,task_directory,restart=True):
        VaspSimulation.__init__(self,task_name,task_directory,restart)

    def config_incar(self):
        VaspSimulation.config_incar(self)

        if self.incar.ibrion not in [1,2]:
            self.incar.ibrion = 2 # set to conjugate gradient method
        self.incar.isif = 2 # relax positions only
        self.incar.write(os.path.join(self.task_directory,'INCAR'))

class VaspCalculatePhonons(VaspSimulation):
    """ calculate the phonons """
    def __init__(self,task_name,task_directory,restart=True):
        VaspSimulation.__init__(self,task_name,task_directory,restart=restart)

    def config_incar(self):
        VaspSimulation.config_incar(self)

        if self.incar.ibrion not in [5,6,7,8]:
            self.incar.ibrion = 6 # frozen phonon mode
            self.incar.potim = 0.015 # displacements
        self.incar.isif = 3
        self.incar.write(os.path.join(self.task_directory,'INCAR'))

class VaspKpointsConvergence(object):
    pass

class VaspEncutConvergence(object):
    """
    Args:
        structure(str or pypospack.crystal.Structure):
        xc(str): exchange correlation functional
        encut_min(int): minimum of the range for the energy cutoff, in eV
        encut_max(int): maximum of the range for the energy cutoff, in eV
        encut_step(int): interval between simulations of energy cutoffs, in eV
        incar_dict(str): parameters for the VASP simulation
        slurm_dict(str): job submission parameters for the VASP simulation
    Attributes:
        poscar(pypospack.io.vasp.Poscar)
        encut_min(int)
        encut_max(int)
        encut_step(int)
        task(dict)
        task_results(dict)
    """
    def __init__(self,directory=None,structure='POSCAR',xc='GGA',
            encut_min=None,encut_max=None,encut_step=25,
            incar_dict=None,slurm_dict=None,full_auto=True):
        
        self.orig_directory = os.getcwd()
        self.filename_results = 'converg.encut.results'
        self.manifest_filename = 'pypospack.manifest.yaml'

        # --- PROCESS ARGUMENTS ---
        # process the dictionary attribute
        if directory is not None:
            if os.path.isabs(directory):
                self.directory = directory
            else:
                self.directory = os.path.join(
                    self.orig_directory,
                    directory)

            if not os.path.exists(self.directory):
                os.makedirs(self.directory)

        # change the location of the result output file
        self.filename_results = os.path.join(
            self.directory,
            self.filename_results)

        self.manifest_filename = os.path.join(
            self.directory,
            self.manifest_filename)
        # check argument 'structure'
        if isinstance(structure,crystal.SimulationCell):
            self.structure_filename = 'POSCAR'
            self.poscar = vasp.Poscar(structure)
            self.poscar.write(structure_filename)
        elif isinstance(structure,str):
            if os.path.isabs(structure):
                self.structure_filename = structure
            else:
                # convert to absolute path
                self.structure_filename = os.path.join(
                    self.orig_directory,
                    structure)
                self.structure_filename = os.path.normpath(
                    self.structure_filename)

            # read the poscar file
            self.poscar = vasp.Poscar()
            self.poscar.read(self.structure_filename)
        else:
            msg_err = "structure must either be an instance of "
            msg_err += "pypospack.crystal.SimulationCell or a string"
            raise ValueError(msg_err)

        # self exchange correlation functional
        self.xc = xc

        # determine energy cutoff
        if any([encut_min is None, encut_max is None]):
            # create a potcar containing the information contained in the 
            # potcar files
            potcar = vasp.Potcar()
            potcar.symbols = self.poscar.symbols
            potcar.xc = xc

            # this is kind of clunky because I have to write the potcar file
            # before i read it.
            # TODO: what I should actually do is read the POTCARS for each
            #       individual symbol and collect the information I need
            potcar.write('POTCAR.tmp')
            potcar.read('POTCAR.tmp')
            os.remove('POTCAR.tmp')

            # set encut min
            if any([isinstance(encut_min,float),
                    isinstance(encut_min,int)]):
                self.encut_min = encut_min
            else:
                self.encut_min = max(potcar.encut_min)

            # set encut max
            if any([isinstance(encut_max,float),
                    isinstance(encut_max,int)]):
                self.encut_max = encut_max
            else:
                self.encut_max = 1.5 * max (potcar.encut_max)
        else:
            self.encut_min = encut_min
            self.encut_max = encut_max

        self.encut_step = encut_step

        # set incut_dict
        if isinstance(incar_dict, dict):
            self.incar_dict = copy.deepcopy(incar_dict)

        # set slurm_dict
        if isinstance(slurm_dict, dict):
            self.slurm_dict = copy.deepcopy(slurm_dict)

        # additional attributes which aren't parameters
        self.tasks = {}
        self.task_results = {}

        if full_auto is True:
            self.do_full_auto()

    def do_full_auto(self):
        if not pathlib.Path(self.manifest_filename).is_file():
            self.create_simulations()
            # write_simulation_manifest()

            self.manifest = SlurmSimulationManifest()
            self.manifest.task_list = copy.deepcopy(self.task_list)
            for task_name in self.task_list:
                task_directory = self.tasks[task_name].task_directory
                #TODO:
                # task_is_created = self.task[task_name].is_created
                # task_is_submitted = self.task[task_name].is_submitted
                # task_is_compled = self.task[task_name].is_completed
                self.manifest.tasks[task_name] = {}
                self.manifest.tasks[task_name]['directory'] = task_directory
                self.manifest.tasks[task_name]['created'] = True
                self.manifest.tasks[task_name]['submitted'] = False
                self.manifest.tasks[task_name]['complete'] = False
            
            self.run_simulations()
            for task_name in self.task_list:
                self.manifest.tasks[task_name]['sim_submitted'] = True
            self.manifest.write(self.manifest_filename)
        else:
            # read_simulation_manifest()
            self.manifest = SlurmSimulationManifest()
            self.manifest.read(self.manifest_filename)
            self.task_list = list(self.manifest.task_list)
            self.tasks = {}

            for task_name in self.task_list:
                task_directory = self.manifest.tasks[task_name]['directory']
                self.tasks[task_name] = VaspSimulation(\
                    task_name=task_name,
                    task_directory=task_directory,
                    restart=True)

                status = self.tasks[task_name].status
                if status == 'CONFIG':
                    pass
                elif status == 'RUN':
                    self.manifest.tasks[task_name]['submitted'] = True
                elif status == 'POST':
                    self.manifest.tasks[task_name]['complete'] = True

            all_sims_complete = all([v['complete'] for k,v in \
                    self.manifest.tasks.items()])

            if all_sims_complete:
                self.post_process()

                # write out everything
                names = ['encut','total_energy','n_atoms','total_energy_per_atom']
                values = []
                for task_name in self.task_list:
                    values.append([self.task_results[task_name][n] for n in names])
                str_out = ','.join(names)+'\n'
                for row in values:
                    str_out += ','.join([str(v) for v in row])+'\n'
                with open(self.filename_results,'w') as f:
                    f.write(str_out)
            else:
                print('simulations not complete')

    def create_simulations(self):
        # local copy
        encut_min = self.encut_min
        encut_max = self.encut_max
        encut_step = self.encut_step

        encuts = np.arange(encut_min,encut_max+1,encut_step).tolist()

        for encut in encuts:
            encut = int(encut)
            str_encut = str(encut)
            incar_dict = copy.deepcopy(self.incar_dict)
            incar_dict['encut'] = encut
            
            task_name = str_encut
            task_directory = os.path.join(self.directory,task_name)
            self.tasks[task_name] = VaspSimulation(
                        task_name=task_name,
                        task_directory=task_directory)
            try:
                self.tasks[task_name].config(
                        poscar=self.structure_filename,
                        incar=incar_dict,
                        xc=self.xc)
            except AttributeError as e:
                s = str(e)
                print('debugging code')
                print(type(self))
                print('\ttask_name:{}'.format(task_name))
                print('\ttask_obj:{}'.format(str(type(self.tasks[task_name]))))
                print('\ttask_obj.potcar:{}'.format(
                        str(type(self.tasks[task_name].potcar))))
                raise

        self.task_list = [str(int(v)) for v in encuts]

    def run_simulations(self):
        for task_name,task in self.tasks.items():
            task.run(job_type='slurm',exec_dict=self.slurm_dict)

    def post_process(self):
        self.task_results = {}
        for task_name,task in self.tasks.items():
            if task.status == 'POST':
                task.postprocess()

                encut = task.outcar.encut
                total_energy = task.outcar.total_energy
                n_atoms = len(task.contcar.atomic_basis)
                total_energy_per_atom = total_energy/n_atoms

                self.task_results[task_name] = {
                    'encut':encut,
                    'total_energy':total_energy,
                    'n_atoms':n_atoms,
                    'total_energy_per_atom':total_energy_per_atom}
         
# ---- temporary classes until i can write a VaspSimulationManager ---

class SimulationManifest(object):
    def __init__(self,filename='pypospack.manifest.yaml'):
        raise NotImplementedError
    def process(self,filename=None):
        raise NotImplementedError
    def read(self,filename=None):
        raise NotImplementedError
    def write(self,filename=None):
        raise NotImplementedError

class SlurmSimulationManifest(SimulationManifest):
    def __init__(self,filename='pypospack.manifest.yaml'):
        self.filename = filename
        self.manifest_dict = {}
        self.task_list = []
        self.tasks = {}

    def read(self,filename=None):
        if filename is not None:
            self.filename = filename

        with open(self.filename,'r') as f:
            manifest = yaml.load(f)

        self.task_list = list(manifest['task_list'])
        self.tasks = dict(manifest['tasks'])

    def write(self,filename):
        manifest = {}
        manifest['task_list'] = self.task_list
        manifest['tasks'] = self.tasks
        if filename is not None:
            self.filename = filename

        with open(self.filename,'w') as f:
            yaml.dump(manifest,f,default_flow_style=False)


