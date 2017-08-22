import os, shutil, subprocess
import ase.build
import pypospack.io.vasp as vasp
import pypospack.crystal as crystal
import pypospack.io.slurm as slurm

#
class WorkflowManager(object):
    def __init__(self):
        pass

    def add_task():

class VaspMinimizeStructure(object):
    def __init__(self):

        self.xc = xc
        self.sim_dir = sim_dir

        self.is_job_pending = True
        self.is_job_starting = True
        self.is_job_submitted = False
        self.is_job_completed = False
        if os.path.exists(self.sim_dir):
            self.is_job_submitted = os.path.exists(self.sim_dir,'job.submitted')
            self.is_job_completed = os.path.exists(self.sim_dir,'job.completed')
            if self.is_job_completed:
                self.postprocesses()
        else:
            os.mkdir(self.sim_dir)
            self.create_simulation()
            self.submit_job()

    def create_simulation(self):
        # initialize input files
        self.vs = vasp.VaspSimulation()

        self.vs.poscar = vasp.Poscar(obj_structure)
        self.vs.incar = vasp.Incar()
        self.vs.potcar = vasp.Potcar()
        self.vs.kpoints = vasp.Kpoints()

        self.vs.xc = self.xc
        self.vs.simulation_directory = self.sim_dir
        self.vs.symbols = self.vs.poscar.symbols

        # configure potcar file
        self.vs.potcar.symbols = self.vs.symbols
        self.vs.potcar.xc = self.vs.xc
        fn_potcar = os.path.join(self.vs.simulation_directory,'POTCAR')
        self.vs.potcar.write(fn_potcar)
        self.vs.potcar.read(fn_potcar)
        self.vs.encut = max(self.vs.potcar.encut_max)
        self.vs.natoms = self.vs.poscar.n_atoms

        # configure incar file
        magmom_0 = 1.0
        self.vs.incar.ismear = 1
        self.vs.incar.sigma = 0.20

        self.vs.incar.ispin = 2
        self.vs.incar.magmom = "{}*{}".format(self.vs.natoms,magmom_0)

        self.vs.incar.ibrion = 2
        self.vs.incar.isif = 3
        self.vs.incar.potim = 0.5
        self.vs.incar.ediffg = -0.001

        self.vs.poscar.write(os.path.join(\
                self.vs.simulation_directory,"POSCAR"))
        self.vs.incar.write(os.path.join(\
                self.vs.simulation_directory,"INCAR"))
        self.vs.kpoints.write(os.path.join(\
                self.vs.simulation_directory,"KPOINTS"))

    def submit_job(self):
        pass

    def postprocess(self):
        pass

class VaspConvergeEncut(object):
    pass

class VaspConvergeKpoints(object):
    pass

class VaspMinimizeStructure(object):
    pass

class VaspCalculateElastic(object):
    pass


def calculate_bulk_properties(sim_dir):
    os.mkdir(sim_dir)

if __name__ == '__main__':

    structures = {}
    structures['Ni_fcc_cubic'] = {'symbols':['Ni'],
                                  'sg':'fcc',
                                  'a0':3.508,
                                  'shape':'cubic'}
    structures['Ni_bcc_cubic'] = {'symbols':['Ni'],
                                  'sg':'bcc', 
                                  'a0':3.508,
                                  'shape':'cubic'}
    structures['Ni_hcp_cubic'] = {'symbols':['Ni'],
                                  'sg':'hcp', 
                                  'a0':3.508,
                                  'shape':'cubic'}
    structures['Ni_dia_cubic'] = {'symbols':['Ni'],
                                  'sg':'diamond', 
                                  'a0':3.508,
                                  'shape':'cubic'}
    structures['Ni_sc_cubic']  = {'symbols':['Ni'],
                                  'sg':'sc',
                                  'a0':3.508,
                                  'shape':'cubic'}
    
    root_dir = os.getcwd()


    for k,v in structures.items():
        sim_dir = os.path.join(root_dir,k)
        
        # add in aliases for space groups here
        if v['sg'] in ['dia']:
            v['sg'] = 'diamond'
        if v['sg'] in ['fcc','bcc','diamond','sc','hcp']:
            if isinstance(v['symbols'],list):
                if len(v['symbols']) == 1:
                    v['symbols'] = v['symbols'][0]
                else:
                    raise KeyError('cannot have more than one symbol in {}'.format(sg))
        else:
            raise KeyError('sg is an unsupported space group, passed in {}'.format(sg))

        obj_poscar = None
        if v['shape'] == 'cubic':
            obj_poscar = vasp.Poscar(\
                    ase.build.bulk(\
                        v['symbols'],
                        v['sg'],
                        a=v['a0'],
                        cubic=True))
        elif v['shape'] == 'ortho':
            obj_poscar= vasp.Poscar(\
                    ase.build.bulk(\
                        v['symbols'],
                        v['sg'],
                        a=v['a0'],
                        ortho=True))
        elif v['shape'] == 'prim':
            obj_poscar = vasp.Poscar(\
                    ase.build.bulk(\
                        v['symbols'].
                        v['sg'],
                        a0=v['a0']))
        else:
            raise KeyError('sg is an unsupported space group, pass in {}'.format(sg))


        VaspCalculateBulkProperties(sim_dir,obj_poscar)

