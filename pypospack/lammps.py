# -*- coding utf-8 -*-
import os, subprocess, shutil, math, time
import numpy as np

import pypospack.io.vasp as vasp
import pypospack.io.crystal as crystal
import pypospack.potential as potential

def make_lammps_structure_file(structure, fname_out, sym_order):
    structure_file = StructureFile()
    structure_file.structure = structure
    structure_file.write(fname_out, sym_order)

class SimulationManager():
    """

    Attributes:
        structure_db (dict)
        variable_names (:obj:list of :obj:str): a list of variable names
        potential (pypospack.potential): an instance of the potential
    """
    def __init__(self):
        self._lmps_sim_obj = {}
        self._structure_db = {}
        self._dir_structure_db = None
        self._dir_lammps_sim_db = None
        self.variable_names = None
        self.variable_dict = None
        self.potential = None

    @property
    def variable_names(self):
        return self._variable_names

    @variable_names.setter
    def variable_names(self, var_names):
        assert type(var_names),list
        self._variable_names = var_names

    @property
    def structure_db(self):
        return self._structure_db

    @structure_db.setter
    def structure_db(self,sdb):
        assert type(sdb), dict
        self._structure_db = sdb

    @property
    def dir_structure_db(self):
        return self._dir_structure_db

    @dir_structure_db.setter
    def dir_structure_db(self, dir_sdb):
        assert type(dir_sdb), str
        self._dir_structure_db = dir_sdb

    @property
    def dir_lmps_sim_db(self):
        return self._dir_lmps_sim_db

    @dir_lmps_sim_db.setter
    def dir_lmps_sim_db(self, dir_lsdb):
        assert type(dir_lsdb),str
        self._dir_lmps_sim_db = dir_lsdb

    @property
    def potential(self):
        return self._potential

    @potential.setter
    def potential(self,pot_obj):
        assert type(pot_obj),isinstance(Potential)
        self._potential = pot_obj

    @property
    def variable_dict(self):
        return self._variable_dict

    def determine_simulations(self):
        for vn in self._variable_names:
            structure = vn.split('.')[0]
            variable  = vn.split('.')[1]

            sim_type = None
            if variable in ['E_min','a0','a1','a2','a3']:
                sim_type = 'E_min_all'
            elif variable in ['c11','c12','c44']:
                sim_type = 'elastic'
            elif variable in ['E_sp','a1_sp','a2_sp','a3_ap']:
                sim_type = 'single_point'
            elif variable in ['E_min_pos','a1_min_pos','a2_min_pos','a3_min_pos']:
                sim_type = 'E_min_pos'
            elif variable in ['n_atoms']:
                sim_type = 'structure_info'
            else:
                raise ValueError("variable[{}] not supported".format(variable))

            sim_name = "{}.{}".format(structure,sim_type)
            if sim_name not in self._lmps_sim_obj.keys():
                self.add_lammps_simulation(sim_name,sim_type,structure)

    def add_lammps_simulation(self,sim_name,sim_type,structure_name):
        sim_dir = sim_name

        if sim_type in ['sp','single_point']:
            self._lmps_sim_obj[sim_name] = SinglePointSimulation(sim_name,sim_dir,structure_name)
        elif sim_type in ['elas', 'elastic']:
            self._lmps_sim_obj[sim_name] = ElasticTensorSimulation(sim_name,sim_dir,structure_name)
        elif sim_type in ['E_min_all']:
            self._lmps_sim_obj[sim_name] = EnergyMinimizationSimulation(sim_name, sim_dir, structure_name)
        elif sim_type in ['structure_info']:
            self._lmps_sim_obj[sim_name] = StructureInfo(sim_name,sim_dir,structure_name)
        elif sim_type in ['E_min_pos']:
            self._lmps_sim_obj[sim_name] = PositionMinimizationSimulation(sim_name, sim_dir, structure_name)
        else:
            raise ValueError('sim_type[{}] is not a supported sim_type')

    def evaluate_parameter_set(self,param_dict):
        assert type(param_dict),dict
        # write potential file
        if type(self._potential) == BuckinghamPotential:
            self._fname_potential = 'potential.mod'
            self._potential.write_potential_file(self._fname_potential,param_dict)
        elif type(self._potential) == EamPotential:
            self._fname_potential = 'potential.mod'
            self._fname_eam = 'eam.alloy'
            self._potential.write_potential_file(self._fname_potential,param_dict)
            self._potential.write_setfl_file()
        else:
            msg_err = "unknown potential[{}]"
            msg_err = msg_err.format(type(self._potential))
            raise ValueError(msg_err)

        # evaluate minimizations first
        for k in self._lmps_sim_obj.keys():
            if k.endswith('E_min_all'):
                s_name      = self._lmps_sim_obj[k].structure_name
                s_fname     = os.path.join(self._dir_structure_db,
                                           self._structure_db[s_name][0])
                ls_template = os.path.join(self._dir_lmps_sim_db,
                                           'E_min_all')
                self._lmps_sim_obj[k].dir_sim_template = ls_template
                self._lmps_sim_obj[k].fname_structure_vasp = s_fname
                self._lmps_sim_obj[k].create()

        # run simulation
        for k in self._lmps_sim_obj.keys():
            if k.endswith('E_min_all'):
                self._lmps_sim_obj[k].run()

        # monitor simulation
        all_sims_done = False
        sims_has_error = False
        while (not all_sims_done) or sims_has_error:
            time.sleep(.05)
            list_is_done    = [self._lmps_sim_obj[k].is_done() 
                                   for k in self._lmps_sim_obj.keys() 
                                   if k.endswith('E_min_all')]
            list_has_errors = [self._lmps_sim_obj[k].has_error() 
                                   for k in self._lmps_sim_obj.keys() 
                                   if k.endswith('E_min_all')]

            if False not in list_is_done:
                all_sims_done = True

            if True in list_has_errors:
                sims_has_error = True

        # post-process
        for k in self._lmps_sim_obj.keys():
            if k.endswith('E_min_all'):
                self._lmps_sim_obj[k].postprocess()

        lattice_parameters = {}
        for k in self._lmps_sim_obj.keys(): 
            if k.endswith('E_min_all'):
                s_name = k.split('.')[0]
                lattice_parameters[s_name]=self._lmps_sim_obj[k].a1

        # evaluate everything else
        for k in self._lmps_sim_obj.keys():
            if not k.endswith('E_min_all'):
                sim_type = k.split('.')[1]
                s_name      = self._lmps_sim_obj[k].structure_name
                s_fname     = os.path.join(self._dir_structure_db,
                                           self._structure_db[s_name][0])
                ls_template = os.path.join(self._dir_lmps_sim_db,
                                           sim_type)
                self._lmps_sim_obj[k].dir_sim_template = ls_template

                # information from structure, this only works because there
                # is only one structure prototype
                # TODO: generalize this code
                self._lmps_sim_obj[k].fname_structure_vasp = s_fname
                self._lmps_sim_obj[k].new_a1 = lattice_parameters['MgO_NaCl']

                # create lammps simulation
                self._lmps_sim_obj[k].create()

        # run simulation
        for k in self._lmps_sim_obj.keys():
            if not k.endswith('E_min_all'):
                self._lmps_sim_obj[k].run()
 
        # monitor simulation
        all_sims_done = False
        sims_has_error = False
        while (not all_sims_done) or sims_has_error:
            time.sleep(.05)
            list_is_done    = [self._lmps_sim_obj[k].is_done() 
                                   for k in self._lmps_sim_obj.keys() 
                                   if not k.endswith('E_min_all')]
            list_has_errors = [self._lmps_sim_obj[k].has_error() 
                                   for k in self._lmps_sim_obj.keys() 
                                   if not k.endswith('E_min_all')]

            if False not in list_is_done:
                all_sims_done = True

            if True in list_has_errors:
                sims_has_error = True

        # post-process
        for k in self._lmps_sim_obj.keys():
            if not k.endswith('E_min_all'):
                self._lmps_sim_obj[k].postprocess()

        self._variable_dict = {}
        for var in self._variable_names:
            v = None
            s_name = var.split('.')[0]
            v_name = var.split('.')[1]
            if v_name in ['E_min','a0','a1','a2','a3']:
                k = "{}.E_min_all".format(s_name)
                if v_name == 'E_min':
                    v = self._lmps_sim_obj[k].total_energy
                elif v_name == 'a0':
                    v = self._lmps_sim_obj[k].a1
                elif v_name == 'a1':
                    v = self._lmps_sim_obj[k].a1
                elif v_name == 'a2':
                    v = self._lmps_sim_obj[k].a2
                elif v_name == 'a3':
                    v = self._lmps_sim_obj[k].a3
            elif v_name in ['c11','c12','c44']:
                k = "{}.elastic".format(s_name)
                if v_name == 'c11':
                    v = self._lmps_sim_obj[k].c11
                elif v_name == 'c12':
                    v = self._lmps_sim_obj[k].c12
                elif v_name == 'c44':
                    v = self._lmps_sim_obj[k].c44
            elif v_name in ['E_sp','a1_sp','a2_sp','a3_sp']:
                k = "{}.single_point".format(s_name)
                if v_name == 'E_sp':
                    v = self._lmps_sim_obj[k].total_energy
                elif v_name == 'a1_sp':
                    v = self._lmps_sim_obj[k].a1
                elif v_name == 'a2_sp':
                    v = self._lmps_sim_obj[k].a2
                elif v_name == 'a3_sp':
                    v = self._lmps_sim_obj[k].a3
            elif v_name in ['E_min_pos','a1_min_pos','a2_min_pos','a3_minPos']:
                k = "{}.E_min_pos".format(s_name)
                if v_name == 'E_min_pos':
                    v = self._lmps_sim_obj[k].total_energy
                elif v_name == 'a1_min_pos':
                    v = self._lmps_sim_obj[k].a1
                elif v_name == 'a2_min_pos':
                    v = self._lmps_sim_obj[k].a2
                elif v_name == 'a3_min_pos':
                    v = self._lmps_sim_obj[k].a3
            elif v_name == 'n_atoms':
                k = "{}.structure_info".format(s_name)
                if v_name == 'n_atoms':
                    v = self._lmps_sim_obj[k].n_atoms
            self._variable_dict[var] = v

class Simulation(base.Simulation):

    def __init__(self, 
                 sim_name,
                 sim_dir,
                 structure_name):
        self._sim_name = sim_name
        self._sim_dir = sim_dir
        self._dir_sim_template = None
        self._structure_name = structure_name
        self._new_a1 = None

        self._fname_structure_vasp = None
        self._fname_lammps_log = 'log.out'
        self._fname_lammps_script = None
        self._fname_potential = 'potential.mod'
        self._fname_run_lammps = 'runsimulation.sh'
        self._fname_lammps_results = 'out.dat'

        self._is_done = False     #boolean flag set if simulation is done
        self._has_error = False   #boolean flag set if simulltion has errors
        self._err_str = None      #string to set if simulation has error

    @property
    def lammps_log_file_filename(self):
        return self._fname_lammps_log

    @lammps_log_file_filename.setter
    def lammps_log_file_filename(self,fname):
        self._fname_lammps_log = fname

    @property
    def structure_name(self):
        return self._structure_name

    @structure_name.setter
    def structure_name(self, name):
        assert type(name), str
        self._structure_name = name

    @property
    def fname_structure_vasp(self):
        return self._fname_structure_vasp

    @fname_structure_vasp.setter
    def fname_structure_vasp(self,fname):
        assert type(fname), str
        self._fname_structure_vasp = fname

    @property
    def dir_sim_template(self):
        return self._dir_sim_template

    @dir_sim_template.setter
    def dir_sim_template(self,dir_t):
        assert type(dir_t),str
        self._dir_sim_template = dir_t

    @property
    def new_a1(self):
        return self._new_a1

    @new_a1.setter
    def new_a1(self,a1):
        assert type(a1), float
        self._new_a1 = a1

    def create(self):
        # remove the directory if it exists
        if os.path.exists(self._sim_dir):
            shutil.rmtree(self._sim_dir)
        
        # copy simulation template
        shutil.copytree(src = self._dir_sim_template,
                        dst = self._sim_dir)

        # create_structure_file
        src = self._fname_structure_vasp
        dst = os.path.join(os.path.curdir,
                           self._sim_dir,
                           'lammps.structure')
        self._create_structure_file(src,dst)
        # copy potential file
        src = self._fname_potential
        dst = os.path.join(os.path.curdir,
                           self._sim_dir,
                           'potential.mod')
        shutil.copyfile(src,dst)

    def _create_structure_file(self,src,dst):
        # read in vasp file
        src_vasp_structure = vasp.Poscar()
        src_vasp_structure.read_file(fname_in=src)

        if self._new_a1 is not None:
            src_vasp_structure.lattice_parameter = self._new_a1

        dst_lmps_structure = StructureFile(obj=src_vasp_structure)
        dst_lmps_structure.write_file(fname_out=dst)

    def run(self):
        self._process = subprocess.Popen('bash ' + self._fname_run_lammps,
                                         shell=True,
                                         cwd=self._sim_dir)
        self._is_done = False
        self._has_error = False

    def has_error(self):
        lines = base.tail(fname = self._fname_lammps_log, 
                          n_lines = 1)
        for line in lines:
            if 'Neighbor list overflow, boost neigh_modify one' in line:
                self._has_error = True
                self._err_str = 'Neighbor list overflow, boost neigh_modify one'
            else:
                self._has_error = False

        return self._has_error

    def is_done(self):
        if self._process.poll() is not None:
            self._is_done = True
        else:
            self._is_done = False
        return self._is_done

    def postprocess(self):
        raise NotImplementedError

class StructureInfo(Simulation):
    def __init__(self,sim_name,sim_dir,structure_name):
        Simulation.__init__(self,sim_name,sim_dir,structure_name)
        self._n_atoms = None
        self._is_done = True
        self._has_error = False

    @property
    def n_atoms(self): return self._n_atoms

    def create(self): pass
    def run(self): pass
    def is_done(self): return True
    def has_error(self): return False
    def postprocess(self):
        src_vasp_structure = vasp.Poscar()
        src_vasp_structure.read_file(self._fname_structure_vasp)
        self._n_atoms = src_vasp_structure.n_atoms

class SinglePointSimulation(Simulation):
    def __init__(self,sim_name,sim_dir,structure_name):
        Simulation.__init__(self,sim_name,sim_dir,structure_name)
        self._a1 = None; self._a2 = None; self._a3 = None
        self._xx = None; self._yy = None; self._zz = None
        self._xy = None; self._xz = None; self._yz = None
        self._pxx = None; self._pyy = None; self._pzz = None
        self._pxy = None; self._pxz = None; self._pyz = None
        self._total_energy  = None

    @property
    def a1(self): return self._a1

    @property
    def a2(self): return self._a2

    @property
    def a3(self): return self._a3

    @property
    def total_energy(self): return self._total_energy

    def postprocess(self):
        fname_results_file = os.path.join(self._sim_dir,self._fname_lammps_results)
        f = open(fname_results_file,'r')
        lines = f.readlines()
        f.close()

        for line in lines:
            # cell simulation variables
            if line.startswith('xx'):
                self._a1 = float(line.split('=')[1])
                self._xx = float(line.split('=')[1])
            elif line.startswith('yy'):
                self._a2 = float(line.split('=')[1])
                self._yy = float(line.split('=')[1])
            elif line.startswith('zz'):
                self._a3 = float(line.split('=')[1])
                self._zz = float(line.split('=')[1])
            elif line.startswith('xy'):
                self._xy = float(line.split('=')[1])
            elif line.startswith('xz'):
                self._xz = float(line.split('=')[1])
            elif line.startswith('yz'):
                self._yz = float(line.split('=')[1])
            # pressure tensor variable
            elif line.startswith('pxx'):
                self._pxx = float(line.split('=')[1])
            elif line.startswith('pyy'):
                self._pyy = float(line.split('=')[1])
            elif line.startswith('pzz'):
                self._pzz = float(line.split('=')[1])
            elif line.startswith('pxy'):
                self._pxy = float(line.split('=')[1])
            elif line.startswith('pxz'):
                self._pxz = float(line.split('=')[1])
            elif line.startswith('pyz'):
                self._pyz = float(line.split('=')[1])
            # energy varaibles
            elif line.startswith('tot_energy'):
                self._total_energy = float(line.split('=')[1])
            elif line.startswith('ecoh'):
                self._energy_per_atom = float(line.split('=')[1])

class EnergyMinimizationSimulation(Simulation):
    def __init__(self,sim_name,sim_dir,structure_name):
        Simulation.__init__(self,sim_name,sim_dir,structure_name)
        self._a1 = None; self._a2 = None; self._a3 = None
        self._xx = None; self._yy = None; self._zz = None
        self._xy = None; self._xz = None; self._yz = None
        self._pxx = None; self._pyy = None; self._pzz = None
        self._pxy = None; self._pxz = None; self._pyz = None
        self._total_energy  = None

    @property
    def a1(self): return self._a1

    @property
    def a2(self): return self._a2

    @property
    def a3(self): return self._a3

    @property
    def total_energy(self): return self._total_energy


    def postprocess(self):
        fname_results_file = os.path.join(self._sim_dir,self._fname_lammps_results)
        f = open(fname_results_file,'r')
        lines = f.readlines()
        f.close()

        for line in lines:
            # cell simulation variables
            if line.startswith('xx'):
                self._a1 = float(line.split('=')[1])
                self._xx = float(line.split('=')[1])
            elif line.startswith('yy'):
                self._a2 = float(line.split('=')[1])
                self._yy = float(line.split('=')[1])
            elif line.startswith('zz'):
                self._a3 = float(line.split('=')[1])
                self._zz = float(line.split('=')[1])
            elif line.startswith('xy'):
                self._xy = float(line.split('=')[1])
            elif line.startswith('xz'):
                self._xz = float(line.split('=')[1])
            elif line.startswith('yz'):
                self._yz = float(line.split('=')[1])
            # pressure tensor variable
            elif line.startswith('pxx'):
                self._pxx = float(line.split('=')[1])
            elif line.startswith('pyy'):
                self._pyy = float(line.split('=')[1])
            elif line.startswith('pzz'):
                self._pzz = float(line.split('=')[1])
            elif line.startswith('pxy'):
                self._pxy = float(line.split('=')[1])
            elif line.startswith('pxz'):
                self._pxz = float(line.split('=')[1])
            elif line.startswith('pyz'):
                self._pyz = float(line.split('=')[1])
            # energy varaibles
            elif line.startswith('tot_energy'):
                self._total_energy = float(line.split('=')[1])
            elif line.startswith('ecoh'):
                self._energy_per_atom = float(line.split('=')[1])

class PositionMinimizationSimulation(Simulation):
    def __init__(self,sim_name,sim_dir,structure_name):
        Simulation.__init__(self,sim_name,sim_dir,structure_name)
        self._a1 = None; self._a2 = None; self._a3 = None
        self._xx = None; self._yy = None; self._zz = None
        self._xy = None; self._xz = None; self._yz = None
        self._pxx = None; self._pyy = None; self._pzz = None
        self._pxy = None; self._pxz = None; self._pyz = None
        self._total_energy  = None

    @property
    def a1(self): return self._a1

    @property
    def a2(self): return self._a2

    @property
    def a3(self): return self._a3

    @property
    def total_energy(self): return self._total_energy


    def postprocess(self):
        fname_results_file = os.path.join(self._sim_dir,self._fname_lammps_results)
        f = open(fname_results_file,'r')
        lines = f.readlines()
        f.close()

        for line in lines:
            # cell simulation variables
            if line.startswith('xx'):
                self._a1 = float(line.split('=')[1])
                self._xx = float(line.split('=')[1])
            elif line.startswith('yy'):
                self._a2 = float(line.split('=')[1])
                self._yy = float(line.split('=')[1])
            elif line.startswith('zz'):
                self._a3 = float(line.split('=')[1])
                self._zz = float(line.split('=')[1])
            elif line.startswith('xy'):
                self._xy = float(line.split('=')[1])
            elif line.startswith('xz'):
                self._xz = float(line.split('=')[1])
            elif line.startswith('yz'):
                self._yz = float(line.split('=')[1])
            # pressure tensor variable
            elif line.startswith('pxx'):
                self._pxx = float(line.split('=')[1])
            elif line.startswith('pyy'):
                self._pyy = float(line.split('=')[1])
            elif line.startswith('pzz'):
                self._pzz = float(line.split('=')[1])
            elif line.startswith('pxy'):
                self._pxy = float(line.split('=')[1])
            elif line.startswith('pxz'):
                self._pxz = float(line.split('=')[1])
            elif line.startswith('pyz'):
                self._pyz = float(line.split('=')[1])
            # energy varaibles
            elif line.startswith('tot_energy'):
                self._total_energy = float(line.split('=')[1])
            elif line.startswith('ecoh'):
                self._energy_per_atom = float(line.split('=')[1])

class ElasticTensorSimulation(Simulation):
    def __init__(self,sim_name,sim_dir,structure_name):
        Simulation.__init__(self,sim_name,sim_dir,structure_name)
        self._c11 = None
        self._c12 = None
        self._c44 = None

    @property
    def c11(self): return self._c11

    @property
    def c12(self): return self._c12

    @property
    def c44(self): return self._c44

    def create(self):
        super(ElasticTensorSimulation,self).create()
        # elastic simulations require some restart variables.
        dst_pot = os.path.join(self._sim_dir,
                               self._fname_potential)
        f = open(dst_pot,'a')
        f.write("# setup minimization style\n")
        f.write("min_style cg\n")
        f.write("min_modify dmax ${dmax} line quadratic\n\n")
        f.write("# setup output\n")
        f.write("thermo 1\n")
        f.write("thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol\n")
        f.write("thermo_modify norm no\n")
        f.close()

    def postprocess(self):
        fname_results_file = os.path.join(self._sim_dir,self._fname_lammps_results)
        f = open(fname_results_file,'r')
        lines = f.readlines()
        f.close()

        for line in lines:
            if line.startswith('c11'):
                self._c11 = float(line.split('=')[1].split()[0])
            elif line.startswith('c12'):
                self._c12 = float(line.split('=')[1].split()[0])
            elif line.startswith('c44'):
                self._c44 = float(line.split('=')[1].split()[0])


class InputFile:
    def __init__(self):
        self.units = "metal"
        self.dimension = 3
        self.boundary = "p p p"
        self.atom_style = "atomic"
        self.atom_modify = "map array"
        self.structure_fname = None

        self.potential_type = "eam/alloy"
        self.eam_alloy_file = "eam.alloy"

    def initializeSimulations(self):
         str_out  = "# written by pyPosMat\n"
         str_out += "# ---- intialize simulation\n"
         str_out += "clear\n"
         str_out += "units {}\n".format(self.units)
         str_out += "dimension {}\n".format(self.dimension)
         str_out += "boundary {}\n".format(self.boundary)
         str_out += "atom_style {}\n".format(self.atom_style)
         str_out += "atom_modify {}\n".format(self.atom_modify)
         return str_out

    def createAtoms(self):
         if self.structure_fname is None:
             err_msg = "LAMMPS structure file not specified."
             raise ValueError(err_msg)

         str_out  = "# ---- create atoms\n"
         str_out += "read_data {}\n".format(self.structure_fname)

    def defineEmpiricalPotential(self):
        str_out  = "# ---- define interatomic potential"
        if self.structure_fname == "eam/alloy":
            str_out += "pair_style eam/alloy\n"
            str_out += "pair_coeff * * {} Ni\n".format(self.eam_alloy_file)
            str_out += "neighbor 2.0 bin\n"
            str_out += "neigh_modify delay 10 check yes\n"
        return str_out

    def defineSettings(self):
        str_out  = "compute eng all pe/atom\n"
        str_out += "compute eatoms all reduce sum c_eng\n"
        return str_out

    def runMinimization(self):
        str_out  = "reset_timestep 0\n"
        str_out += "fix 1 all box/relax iso 0.0 vmax 0.001\n"
        str_out += "thermo 10\n"
        str_out += "thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms\n"
        str_out += "min_style cg\n"
        str_out += "minimize 1e-25 1e-25 5000 10000\n"
        return str_out

    def defineVariables(self):
        str_out  = "# ---- define variables"
        str_out += "variable natoms equal \"count(all)\"\n"
        str_out += "variable tot_energy equal \"c_eatoms\"\n"
        str_out += "variable length_x equal \"lx/4\"\n"
        str_out += "variable length_y equal \"ly/4\"\n"
        str_out += "variable length_z equal \"lz/4\"\n"
        str_out += "variable ecoh equal \"pe/atoms\"\n"
        return str_out

    def output(self):
        str_out  = "# ---- output\n"
        str_out += "print \"pyPosMat output section\"\n"
        str_out += "print \"tot_energy = ${tot_energy}\"\n"
        str_out += "print \"num_atoms = ${natoms}\"\n"
        str_out += "print \"latt_const_a = ${length_x}\"\n"
        str_out += "print \"latt_const_b = ${length_y}\"\n"
        str_out += "print \"latt_const_c = ${length_z}\"\n"
        str_out += "print \"ecoh = ${ecoh}\"\n"
        str_out += "print \"lammps_sim_done\"\n"
        return str_out
  
class StructureFile(base.Structure):
    def __init__(self,obj=None):
        base.Structure.__init__(self,obj)
        
    def write_file(self, fname_out, symbol_list=None, atom_style=None):
        if symbol_list is None:
            symbol_list = self.symbols

        if atom_style is None:
            atom_style = 'charge'

        total_number_of_atoms      = self.n_atoms
        total_number_of_atom_types = len(self.symbols)
        a0 = self.lattice_parameter

        xlo                        = 0.0
        xhi                        = self.h_matrix[0,0] * a0
        ylo                        = 0.0
        yhi                        = self.h_matrix[1,1] * a0
        zlo                        = 0.0
        zhi                        = self.h_matrix[2,2] * a0
        xy                         = self.h_matrix[0,1] * a0
        xz                         = self.h_matrix[0,2] * a0
        yz                         = self.h_matrix[1,2] * a0

        file = open(fname_out,'w')
        file.write("# {}\n".format(symbol_list))
        file.write("\n")
        file.write("{} atoms\n".format(total_number_of_atoms))
        file.write("{} atom types\n".format(total_number_of_atom_types))
        file.write("\n")
        file.write("{:10.4f} {:10.4f} xlo xhi\n".format(xlo, xhi))
        file.write("{:10.4f} {:10.4f} ylo yhi\n".format(ylo, yhi))
        file.write("{:10.4f} {:10.4f} zlo zhi\n".format(zlo, zhi))
        file.write("\n")
        file.write("{:10.4f} {:10.4f} {:10.4f} xy xz yz\n".format(xy,xz,yz))
        file.write("\n")
        file.write("Atoms\n")
        file.write("\n")

        atom_id = 1
        for i_symbol, symbol in enumerate(symbol_list):
            for i_atom, atom in enumerate(self.atoms):
                if (atom.symbol == symbol):
                    chrg = 1.  # dummy variable
                    posx = self.h_matrix[0,0]*atom.position[0]*a0
                    posy = self.h_matrix[1,1]*atom.position[1]*a0
                    posz = self.h_matrix[2,2]*atom.position[2]*a0
                    if atom_style == 'atomic':
                        str_out = "{} {} {:10.4f} {:10.4f} {:10.4f}\n"
                        str_out = str_out.format(atom_id, 
                                                 i_symbol + 1, 
                                                 posx,posy,posz)
                    elif atom_style == 'charge':
                        str_out = "{} {} {:10.4f} {:10.4f} {:10.4f} {:10.4f}\n"
                        str_out = str_out.format(atom_id,
                                                 i_symbol + 1,
                                                 chrg, posx,posy,posz)
                    file.write(str_out)
                    atom_id += 1
        file.close()

# various helper functions

def checkLammpsSimulationDone(fname_lmps_log):
  lines = pyflamestk.base.tail(fname = fname_lmps_log, n_lines = 2)
  if 'lammps_sim_done' in  lines:
    return True
  else:
    return False

def checkAllLammpsSimulationsDone(sim_directories):
  simIsDone = [False for i in sim_directories]
  for idx, dir in enumerate(sim_directories):
    fname = "{}/out.dat".format(dir)
    simIsDone[idx] = checkLammpsSimulationDone(fname)
  if False in simIsDone:
    return False
  else:
    return True

def checkLammpsSimulationForError(fname_lmps_log):
  lines = pyflamestk.base.tail(fname = fname_lmps_log, n_lines = 1)
  for line in lines:
    if 'Neighbor list overflow, boost neigh_modify one' in line:
      return True
  else:
    return False

def checkAllLammpsSimulationsForErrors(sim_directories):
  simHasError = False
  for dir in sim_directories:
    fname = "{}/out.dat".format(dir)
    simHasError = checkLammpsSimulationForError(fname)
    if (simHasError):
      return True
  return False

def createLammpsSimulation(sim_name,
                           fname_structure,
                           fname_sim_template,
                           fname_potential = None,
                           fname_eam = None):
  dir_name = sim_name

  if os.path.exists(dir_name):
    shutil.rmtree(dir_name)

  shutil.copytree(fname_sim_template,
                  "{}".format(dir_name))

  shutil.copyfile(fname_structure,
                  "{}/lammps.structure".format(dir_name))

  if not(fname_potential == None):
    pass

  if not(fname_eam == None):
    shutil.copyfile(fname_eam,
                   "{}/eam.alloy".format(dir_name))

  for dir in sim_directories:
    p = subprocess.Popen("./runsimulation.sh", shell=True, cwd=dir)

  # check if simulations are done
  sims_finished = False
  while not(sims_finished):
    sims_finished = checkAllLammpsSimulationsDone(sim_directories = sim_directories)
    time.sleep(0.05)

def getCohesiveEnergy(fname_lammps_out):
  ecoh=0
  f = open(fname_lammps_out,'r')
  lines = f.readlines()
  for line in lines:
    if line.startswith("ecoh ="):
      ecoh = float(line.split('=')[1].strip())
  f.close()
  return ecoh

def getElasticComponent(fname_lammps_out, component):
  retval = 0
  f = open(fname_lammps_out, 'r')
  lines = f.readlines()
  for line in lines:
    if line.startswith("{} =".format(component)):
      retval = line.split('=')[1].strip().split(' ')[0]
      retval = float(retval)
  f.close()
  return retval

def getLatticeParameter(fname_lammps_out):
  alat = 0
  f = open(fname_lammps_out, 'r')
  lines = f.readlines()
  for line in lines:
    if line.startswith("latt_const_a"):
      alat = float(line.split('=')[1].strip())
  f.close()
  return alat

def getPressure(fname_lammps_out, type = 'total'):
  f = open(fname_lammps_out, 'r')
  lines = f.readlines()
  pressure_xx = 0
  pressure_yy = 0
  pressure_zz = 0
  pressure_total = 0
  for line in lines:
    if line.startswith("pressure_total"):
      pressure_total = float(line.split('=')[1].strip())
    elif line.startswith("pressure_xx"):
      pressure_xx = float(line.split('=')[1].strip())
    elif line.startswith("pressure_yy"):
      pressure_yy = float(line.split('=')[1].strip())
    elif line.startswith("pressure_zz"):
      pressure_zz = float(line.split('=')[1].strip())
    else:
      pass
  if   type == 'total':
    return pressure_total
  elif type == 'xx':
    return pressure_xx
  elif type == 'yy':
    return pressure_yy
  elif type == 'zz':
    return pressure_zz
  elif type == 'all':
    return [pressure_xx, pressure_yy, pressure_zz, pressure_total]
  else:
    print("pyflamestk.lammps -> do not understand pressure_type: {}".format(type))

