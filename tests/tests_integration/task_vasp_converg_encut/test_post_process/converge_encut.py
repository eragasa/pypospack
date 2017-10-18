import yaml, copy, pathlib
import numpy as np
import pypospack.io.vasp as vasp
import pypospack.crystal as crystal
import pypospack.task.vasp as tsk_vasp
import os

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

    def write(self,filename):
        manifest = {}
        manifest['task_list'] = self.task_list
        manifest['tasks'] = self.tasks
        if filename is not None:
            self.filename = filename

        with open(self.filename,'w') as f:
            yaml.dump(manifest,f,default_flow_style=False)

class VaspConvergenceEnergyCutoff(object):
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
    def __init__(self,structure='POSCAR',xc='GGA',
            encut_min=None,encut_max=None,encut_step=25,
            incar_dict=None,slurm_dict=None,full_auto=True):

        # check argument 'structure'
        if isinstance(structure,crystal.SimulationCell):
            self.structure_filename = 'POSCAR'
            self.poscar = vasp.Poscar(structure)
            self.poscar.write(structure_filename)
        elif isinstance(structure,str):
            self.structure_filename = structure
            self.poscar = vasp.Poscar()
            self.poscar.read(structure)
        else:
            msg_err = "structure must either be an instance of "
            msg_err += "pypospack.crystal.SimulationCell or a string"
            raise ValueError(msg_err)

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
        if not pathlib.Path('pypospack.manifest.yaml').is_file():
            self.create_simulations()
            self.manifest = SlurmSimulationManifest()
            self.manifest.task_list = copy.deepcopy(self.task_list)
            for task_name in self.task_list:
                self.manifest.tasks[task_name] = {}
                self.manifest.tasks[task_name]['created'] = True
                self.manifest.tasks[task_name]['submitted'] = False
                self.manifest.tasks[task_name]['complete'] = False
            
            self.run_simulations()
            for task_name in self.task_list:
                self.manifest.tasks[task_name]['sim_submitted'] = True

            self.manifest.write('pypospack.manifest.yaml')
        else:
            self.manifest = SlurmSimulationManifest()
            self.manifest.read('pypospack.manifest.yaml')

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
            self.tasks[str_encut] = tsk_vasp.VaspSimulation(
                        task_name = str_encut,
                        task_directory = str_encut)
            self.tasks[str_encut].config(
                    poscar=self.structure_filename,
                    incar=incar_dict,
                    xc=self.xc)

        self.task_list = [str(int(v)) for v in encuts]

    def run_simulations(self):
        for task_name,task in self.tasks.items():
            task.run(job_type='slurm',exec_dict=self.slurm_dict)

if __name__ == '__main__':
    structure_filename = 'MgO_NaCl_prim.vasp'
    xc = 'GGA'
    full_auto=True
    incar_dict = {}
    incar_dict['ismear'] = 0
    incar_dict['sigma'] = 0.05
    incar_dict['ispin'] = 1

    slurm_dict = {}
    slurm_dict['email'] = 'eragasa@ufl.edu'
    slurm_dict['qos'] = 'phillpot'
    slurm_dict['ntasks'] = 16
    slurm_dict['time'] = "1:00:00"

    encut_conv = tsk_vasp.VaspEncutConvergence(
            structure=structure_filename,xc='GGA',
            incar_dict=incar_dict,slurm_dict=slurm_dict,
            full_auto=full_auto)
    exit()
    for encut in np.arange(encut_min,encut_max+1,encut_step):
        encut = int(encut)
        str_encut = str(encut)
        tasks[str_encut] = tsk_vasp.VaspSimulation(
		    task_name = str_encut,
		    task_directory = str_encut,
		    restart=is_restart)

        incar_config['encut']=encut
        if tasks[str_encut].status == 'INIT':
            tasks[str_encut].config(poscar=structure_filename,incar=incar_config,xc=xc)

        if tasks[str_encut].status == 'CONFIG':
            tasks[str_encut].run(\
                job_type='slurm',
                exec_dict=slurm_dict)

        if tasks[str_encut].status == 'RUN':
            pass

        if tasks[str_encut].status == 'POST':
            tasks[str_encut].postprocess()
            total_energy = tasks[str_encut].outcar.total_energy
            n_atoms = len(tasks[str_encut].contcar.atomic_basis)
            total_energy_per_atom = total_energy / n_atoms

            task_results[str_encut] = {}
            task_results[str_encut]['encut'] = tasks[str_encut].outcar.encut
            task_results[str_encut]['total_energy'] = total_energy
            task_results[str_encut]['n_atoms'] = n_atoms
            task_results[str_encut]['total_energy_per_atom'] = total_energy / n_atoms

    for i,encut in enumerate(np.arange(encut_min,encut_max+1,encut_step)):
        encut = int(encut)
        str_encut = str(encut)
        str_old_encut = str(int(encut-25))
        if i == 0:
            print(str_encut,
                task_results[str_encut]['encut'],
                task_results[str_encut]['total_energy_per_atom']) 
        else:
            print(str_encut,
                task_results[str_encut]['encut'],
                task_results[str_encut]['total_energy_per_atom'],
                100*(task_results[str_encut]['total_energy_per_atom']-task_results[str_old_encut]['total_energy_per_atom']))
         

