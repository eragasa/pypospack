# -*- coding: utf-8 -*- 
""" Implementation of LammpsSimulation abstract and implemented classes 

This module implements simulation tasks.  An abstract class is implemented in
LAMMPS simulation tasks, this class should be subclassed for new implementations
requiring LAMMPS simulations.  Tasks which do not require a LAMMPS simulation 
should subclass the pypospack.io.base.Task instead.

Attributes:
    atom_style_list(:obj:`list` of :obj:`str`)
    potential_map ist: module level variable to indicate atom styles for LAMMPS


Todo:
    * nothing as this point
"""
import os, copy, importlib, subprocess
import numpy as np
import pypospack.potential as potential
import pypospack.io.vasp as vasp
import pypospack.io.lammps as lammps
from pypospack.task import Task

atom_style_list = ['charge','atomic']

potential_map = {\
        'buckingham':['pypospack.potential','Buckingham'],
        'eam':['pypospack.potential','EmbeddedAtomModel'],
        'tersoff':['pypospack.potential','Tersoff']}

class LammpsSimulationError(Exception):
    """Error class for dealing with LAMMPS simulation issues"""
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)

class LammpsSimulation(Task):
    """ Calculates cohesive energy

    This is an abstract data class which defines the attributes and methods 
    necessary to interact with the Workflow manager.  A default implementation
    has been created for this class

    Args:
        task_name(str): unique id for the task name being defined
        task_directory(str): the directory where this task will create
            input and output files for LAMMPS.

    Attributes:
        potential(pypospack.potential.Potential): the potential class
        structure(pypospack.io.vasp.Poscar): the structure
        config(:obj:'list' of :obj:'str'): a list of attributes required to 
            configure this LAMMPS task.
        config_map(dict):
        potential_map(dict):
        results(dict): results of the simulation

    """
    def __init__(self,task_name,task_directory):
        Task.__init__(self,task_name,task_directory)
        if self.status == 'INIT':
            self.potential = None
            self.structure = None
            self.results = None
            additional_config_dict = {\
                'structure':self.set_structure,
                'potential':self.set_potential
                }
            additional_ready_dict = {}
            additional_run_dict = {}
            additional_post_dict = {}
            
            self.config_dict.update(additional_config_dict)
            self.ready_dict.update(additional_config_dict)
            self.run_dict.update(additional_config_dict)
            self.post_dict.update(additional_post_dict)

    def restart(self):
        raise NotImplementedError

    def send_config_info(self,config_dict):
        for k,v in config_dict.items():
            self.config_dict[k](v)

    def config(self, structure=None,potential=None):
        """ configure this class

        Args:
            structure (dict): a dictionary of structural information
            potential (dict): a dictionary of potential information
        """
        if structure is not None:
            self.structure_dict = structure
        if potential is not None:
            self.potential_dict = potential

        self.structure_name = self.structure_dict['name']
        self.structure = vasp.Poscar()
        self.structure.read(self.structure_dict['filename'])
        self.symbols = self.structure.symbols

        pot_type = self.potential_dict['potential_type']
        pot_syms = self.potential_dict['symbols']
        pot_param = self.potential_dict['params']
        self.configure_potential(\
                pot_type = pot_type,
                symbols = pot_syms)
        self.param_dict = copy.deepcopy(pot_param)
        self.status = 'CONFIG'

    def ready(self):
        self.status = 'READY'

    def run(self, is_mpi=False):
        self.status = 'RUN'
        self.write_lammps_input_file()
        self.write_potential_files(param_dict=self.param_dict)
        self.write_structure_file()

        # choose lammps binary and get it from environment variable
        if is_mpi:
            lmps_bin = os.environ.get('LAMMPS_MPI_BIN')
        elif not is_mpi:
            lmps_bin = os.environ.get('LAMMPS_SERIAL_BIN')

        cmd_str = '{} -i lammps.in > lammps.out'.format(lmps_bin)

        # change context directory
        orig_dir = os.getcwd()
        os.chdir(self.task_directory)
        try:
            subprocess.call(cmd_str,shell=True)
        finally:
            os.chdir(orig_dir)

        self.status = 'POST'

    def write_lammps_input_file(self,filename='lammps.in'):
        """ writes LAMMPS input file 
        
        Args:
            filename (str): name of the input file for LAMMPS. Default is
                'lammps.in'.
        """
        str_out = self.lammps_input_file_to_string()
        filename = os.path.join(self.task_directory,filename)
        with open(filename,'w') as f:
            try:
                f.write(str_out)
            except:
                print('str_out:{}'.format(str_out))
                raise

    def lammps_input_file_to_string(self):
        """ string for the LAMMPS input file """

        str_out = "".join([\
                self._lammps_input_initialization_section(),
                self._lammps_input_create_atoms(),
                self._lammps_input_define_potential(),
                self._lammps_input_run_minimization(),
                self._lammps_input_out_section()])
        return(str_out)

    def write_potential_files(self,param_dict=None,filename='potential.mod'):
        filename = os.path.join(self.task_directory,filename)
    
        if param_dict is not None:
            try:
                str_out = self.potential.to_string(param_dict)
            except KeyError as ke:
                s = str(ke)
                print("keyerror:{}".format(s))
                print("param_dict:",param_dict)
                print("self.param_dict:",self.param_dict)
                raise

            with open(filename,'w') as f:
                f.write(str_out)
        else:
            for src_filename in self.potential_files:
                shutil.copyfile(src=self.potential_file,
                                dst=filename)
            
    def write_structure_file(self,filename='structure.dat'):
        filename = os.path.join(self.task_directory,filename)

        if self.potential.is_charge:
            atom_style = 'charge'
        else:
            atom_style = 'atomic'
        symbol_list = self.symbols
    
        # instatiate using lammpsstructure file
        lammps_structure = lammps.LammpsStructure(\
                obj=self.structure)
        lammps_structure.write(\
                filename=filename,
                symbol_list=symbol_list,
                atom_style=atom_style)

    def postprocess(self):
        filename = os.path.join(self.task_directory,'lammps.out')
        with open(filename,'r') as f:
            lines = f.readlines()

        lammps_result_names = ['tot_energy','num_atoms',
                        'xx','yy','zz','xy','xz','yz',
                        'tot_press','pxx','pyy','pzz','pxy','pxz','pyz']

        self.results = {}
        for i,line in enumerate(lines):
            for name in lammps_result_names:
                if '{} = '.format(name) in line:
                    try:
                        self.results[name] = \
                                float(line.split('=')[1].strip())
                    except:
                        print('name:{}'.format(name))
                        print('line:{}'.format(line.strip()))
                        raise

        try:
            self.results['ecoh'] = self.results['tot_energy'] / self.results['num_atoms']
        except KeyError as e:
            print(e)
    def req_config(self):
        """ this should method should be overridden, not inherited """
        return list(self.config_dict.keys())

    def send_config(self,config_dict):
        """ this should method should be overridden, not inherited """
        for k,v in config_dict.items():
            self.config_dict[k](v)
        self.configure()
        self.status = 'CONFIG'

    def req_ready(self):
        """ this should method should be overridden, not inherited """
        return list(self.ready_dict.keys())

    def send_ready(self,ready_dict):
        """ this should method should be overridden, not inherited """
        for k,v in ready_dict:
            try:
                self.ready_dict[k](v)
            except:
                raise

        if None not in [v for k,v in self.ready_dict]:
            self.status = 'READY'

    def req_run(self):
        """ this should method should be overridden, not inherited """
        return list(self.run_dict.keys())

    def send_run(self,run_dict):
        """ this should method should be overridden, not inherited """
        self.status = "RUN"

    def req_post(self):
        """ this should method should be overridden, not inherited """
        return list(self.post_dict.keys())

    def send_post(self,post_dict):
        """ this should method should be overridden, not inherited """
        self.status = "POST"

    def req_done(self):
        """ this should method should be overridden, not inherited """
        return list(self.done_dict.keys())

    def send_done(self,done_dict):
        """ this should method should be overridden, not inherited """
        self.status = "DONE"

    # below here are effectively helper functions
    def set_structure(self,structure_dict):
        self.structure_name = structure_dict['name']
        self.structure_filename = structure_dict['filename']

    def set_potential(self,potential_dict):
        self.potential_type = potential_dict['potential_type']
        self.potential_symbols = potential_dict['symbols']
        self.param_dict = copy.deepcopy(potential_dict['params'])

    def configure_potential(self,pot_type,symbols):
        module_name = potential_map[pot_type][0]
        class_name = potential_map[pot_type][1]

        try:
            module = importlib.import_module(module_name)
            klass = getattr(module,class_name)
            self.potential = klass(symbols)
        except:
            raise

    # private functions for building lammps input files
    def _lammps_input_initialization_section(self):
        if self.potential.is_charge:
            atom_style = 'charge'
        else:
            atom_style = 'atomic'
        str_out = (
            '# ---- initialize simulations\n'
            'clear\n'
            'units metal\n'
            'dimension 3\n'
            'boundary p p p\n'
            'atom_style {}\n'
            'atom_modify map array\n'
              ).format(atom_style)
        return str_out

    def _lammps_input_create_atoms(self):
        str_out = (
            '# ---- create atoms\n'
            'read_data structure.dat\n')

        return str_out

    def _lammps_input_define_potential(self):
        str_out = (
            '# ---- define interatomic potential\n'
            'include potential.mod\n')
        return str_out

    def _lammps_input_run_minimization(self):
        str_out = (
            '# ---- define settings\n'
            'compute eng all pe/atom\n'
            'compute eatoms all reduce sum c_eng\n'
            '# ---- run minimization\n'            
            'reset_timestep 0\n'
            'fix 1 all box/relax iso 0.0 vmax 0.001\n'
            'thermo 10\n'
            'thermo_style custom step pe lx ly lz xy xz yz press pxx pyy pzz pxy pxz pyz c_eatoms\n'
            # 'thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms\n'
            'min_style cg\n'
            'minimize 1e-25 1e-25 5000 10000\n'
            )
        return str_out

    def _lammps_input_out_section(self):
        str_out = (
            '# ---- define output variables ----\n'
            'variable natoms equal "count(all)"\n'
            'variable tot_energy equal "c_eatoms"\n'
            'variable length_x equal "lx"\n'
            'variable length_y equal "ly"\n'
            'variable length_z equal "lz"\n'
            'variable tilt_xy equal "xy"\n'
            'variable tilt_xz equal "xz"\n'
            'variable tilt_yz equal "yz"\n'
            'variable tot_press equal "press"\n'
            'variable press_xx equal "pxx"\n'
            'variable press_yy equal "pyy"\n'
            'variable press_zz equal "pzz"\n'
            'variable press_xy equal "pxy"\n'
            'variable press_xz equal "pxz"\n'
            'variable press_yz equal "pyz"\n'
            '\n'
            '# ---- output ----\n'
            'print \"pypospack:output_section:begin\"\n'
            'print \"tot_energy = ${tot_energy}\"\n'
            'print \"num_atoms = ${natoms}"\n'
            'print \"xx = ${length_x}\"\n'
            'print \"yy = ${length_y}\"\n'
            'print \"zz = ${length_z}\"\n'
            'print \"xy = ${tilt_xy}\"\n'
            'print \"xz = ${tilt_xz}\"\n'
            'print \"yz = ${tilt_yz}\"\n'
            'print \"tot_press = ${tot_press}\"\n'
            'print \"pxx = ${press_xx}\"\n'
            'print \"pyy = ${press_yy}\"\n'
            'print \"pzz = ${press_zz}\"\n'
            'print \"pxy = ${press_xy}\"\n'
            'print \"pxz = ${press_xz}\"\n'
            'print \"pyz = ${press_yz}\"\n'
            'print \"pypospack:output_section:done\"\n'
            'print \"pypospack:lammps_sim:done\"\n'
                  )
        return str_out

class LammpsStructuralMinimization(LammpsSimulation):
    """ Class for LAMMPS structural minimization

    This data class defines additional attributes and methods necessary to 
    interact with the Workflow manager.

    Args:
        task_name(str): unique id for the task name being define
        task_directory(str): the directory where this task will create
            input and output files for LAMMPS

    Attributes:
        config
        config_map
    """
    def __init__(self,task_name,task_directory):
        LammpsSimulation.__init__(self,task_name,task_directory)

    def postprocess(self):
        LammpsSimulation.postprocess(self)

    def lammps_input_file_to_string(self):
        str_out = "".join([\
                self._lammps_input_initialization_section(),
                self._lammps_input_create_atoms(),
                self._lammps_input_define_potential(),
                self._lammps_input_run_minimization(),
                self._lammps_input_out_section()])
        return(str_out)

    def _lammps_input_run_minimization(self):
        str_out = (
            '# ---- define settings\n'
            'compute eng all pe/atom\n'
            'compute eatoms all reduce sum c_eng\n'
            '# ---- run minimization\n'            
            'reset_timestep 0\n'
            'fix 1 all box/relax iso 0.0 vmax 0.001\n'
            'thermo 10\n'
            'thermo_style custom step pe lx ly lz xy xz yz press pxx pyy pzz pxy pxz pyz c_eatoms\n'
            # 'thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms\n'
            'min_style cg\n'
            'minimize 1e-25 1e-25 5000 10000\n'
            )
        return str_out

class LammpsPositionMinimization(LammpsSimulation):
    """ Class for LAMMPS position minimization

    This class sets up, runs and processes the relaxation of atomic positions
    to find the lowest energy structure in the local basin of attraction.  The
    simulation cell is frozen.

    Args:
        task_name (str): the name of the task
        task_directory (str): the directory of the task
    """
    def __init__(self,task_name,task_directory):
        LammpsSimulation.__init__(self,task_name,task_directory)

    def lammps_input_file_to_string(self):
        str_out = "".join([\
                self._lammps_input_initialization_section(),
                self._lammps_input_create_atoms(),
                self._lammps_input_define_potential(),
                self._lammps_input_position_minimization(),
                self._lammps_input_out_section()])
        return(str_out)

    def _lammps_input_position_minimization(self):
        str_out = (\
            "# ---- define settings\n"
            "compute eng all pe/atom\n"
            "compute eatoms all reduce sum c_eng\n"
            "# ---- run minimization\n"
            "reset_timestep 0\n"
            "thermo 1\n"
            "thermo_style custom step pe lx ly lz xy xz yz press pxx pyy pzz pxy pxz pyz c_eatoms\n"
            "min_style cg\n"
            "minimize 1e-20 1e-20 1000 100000\n")
        return str_out

class LammpsElasticCalculation(LammpsSimulation):
    """ Class for LAMMPS elastic calculation

    This class sets up, runs and processes the relaxation of atomic poistions
    to find the lowest energy struction in the local basin of attraction. The
    simulation cell is then deformed and atomic positions relaxed to calculate
    the energy differences to calculate elements of the elastic tensor

    This class is based on the LAMMPS script written by Aiden Thompson

    Args:
        task_name (str): the name of the task
        task_directory (str): the directory of the task
    """
    def __init__(self,task_name,task_directory):
        LammpsSimulation.__init__(self,task_name,task_directory)

    def write_potential_files(self,param_dict=None,filename='potential.mod'):
        LammpsSimulation.write_potential_files(\
                self,
                param_dict,
                filename=filename)

        str_out = (\
                "# setup minimization style\n"
                "min_style cg\n"
                "min_modify dmax ${dmax} line quadratic\n"
                "# setup output\n"
                "thermo 1\n"
                "thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol\n"
                "thermo_modify norm no\n")
        filename = os.path.join(self.task_directory,filename)
        with open(filename,'a') as f:
            f.write(str_out)

    def write_lammps_input_file(self,filename='lammps.in'):
        """ writes LAMMPS input file

        This method is modified from the LammpsSimulation template due to 
        the need for the multiple files.
        
        Args:
            filename (str): name of the input file for LAMMPS. Default is
                'lammps.in'.
        Attributes:
        
        """

        str_out = self.lammps_input_file_to_string()
        filename = os.path.join(self.task_directory,filename)
        with open(filename,'w') as f:
            f.write(str_out)

        str_out = self.lammps_init_mod_to_string()
        filename = os.path.join(self.task_directory,'init.mod')
        with open(filename,'w') as f:
            f.write(str_out)

        str_out = self.lammps_displace_mod_to_string()
        filename = os.path.join(self.task_directory,'displace.mod')
        with open(filename,'w') as f:
            f.write(str_out)

    def lammps_input_file_to_string(self):
        str_out = (
            "include init.mod\n"
            "include potential.mod\n"
            "# ---- Compute initial state\n"
            "fix 3 all box/relax aniso 0.0\n"
            "minimize ${etol} ${ftol} ${maxiter} ${maxeval}\n"
            "\n"
            "variable tmp equal pxx\n"
            "variable pxx0 equal ${tmp}\n"
            "variable tmp equal pyy\n"
            "variable pyy0 equal ${tmp}\n"
            "variable tmp equal pzz\n"
            "variable pzz0 equal ${tmp}\n"
            "variable tmp equal pyz\n"
            "variable pyz0 equal ${tmp}\n"
            "variable tmp equal pxz\n"
            "variable pxz0 equal ${tmp}\n"
            "variable tmp equal pxy\n"
            "variable pxy0 equal ${tmp}\n"
            "\n"
            "variable tmp equal lx\n"
            "variable lx0 equal ${tmp}\n"
            "variable tmp equal ly\n"
            "variable ly0 equal ${tmp}\n"
            "variable tmp equal lz\n"
            "variable lz0 equal ${tmp}\n"
            "\n"
            " # ---- define the derivatives w.r.t. strain components\n"
            "variable d1 equal -(v_pxx1-${pxx0})/(v_delta/v_len0)*${cfac}\n"
            "variable d2 equal -(v_pyy1-${pyy0})/(v_delta/v_len0)*${cfac}\n"
            "variable d3 equal -(v_pzz1-${pzz0})/(v_delta/v_len0)*${cfac}\n"
            "variable d4 equal -(v_pyz1-${pyz0})/(v_delta/v_len0)*${cfac}\n"
            "variable d5 equal -(v_pxz1-${pxz0})/(v_delta/v_len0)*${cfac}\n"
            "variable d6 equal -(v_pxy1-${pxy0})/(v_delta/v_len0)*${cfac}\n"
            "\n"
            "# ---- write restart files\n"
            "unfix 3\n"
            "write_restart restart.equil\n"
            "# ---- uxx Perturbation\n"
            "variable dir equal 1\n"
            "include displace.mod\n"
            "# ---- uyy Perturbation\n"
            "variable dir equal 2\n"
            "include displace.mod\n"
            "# ---- uzz Perturbation\n"
            "variable dir equal 3\n"
            "include displace.mod\n"
            "# ---- uyz Perturbation\n"
            "variable dir equal 4\n"
            "include displace.mod\n"
            "# ---- uxz Perturbation\n"
            "variable dir equal 5\n"
            "include displace.mod\n"
            "# ---- uxy Perturbation\n"
            "variable dir equal 6\n"
            "include displace.mod\n"
            "\n"
            "# ---- Output final values\n"
            "variable C11all equal ${C11}\n"
            "variable C22all equal ${C22}\n"
            "variable C33all equal ${C33}\n"
            "variable C12all equal 0.5*(${C12}+${C21})\n"
            "variable C13all equal 0.5*(${C13}+${C31})\n"
            "variable C23all equal 0.5*(${C23}+${C32})\n"
            "variable C44all equal ${C44}\n"
            "variable C55all equal ${C55}\n"
            "variable C66all equal ${C66}\n"
            "variable C14all equal 0.5*(${C14}+${C41})\n"
            "variable C15all equal 0.5*(${C15}+${C51})\n"
            "variable C16all equal 0.5*(${C16}+${C61})\n"
            "variable C24all equal 0.5*(${C24}+${C42})\n"
            "variable C25all equal 0.5*(${C25}+${C52})\n"
            "variable C26all equal 0.5*(${C26}+${C62})\n"
            "variable C34all equal 0.5*(${C34}+${C43})\n"
            "variable C35all equal 0.5*(${C35}+${C53})\n"
            "variable C36all equal 0.5*(${C36}+${C63})\n"
            "variable C45all equal 0.5*(${C45}+${C54})\n"
            "variable C46all equal 0.5*(${C46}+${C64})\n"
            "variable C56all equal 0.5*(${C56}+${C65})\n"
            "\n"
            "print \"c11 = ${C11all} ${cunits}\"\n"
            "print \"c22 = ${C22all} ${cunits}\"\n"
            "print \"c33 = ${C33all} ${cunits}\"\n"
            "print \"c12 = ${C12all} ${cunits}\"\n"
            "print \"c13 = ${C13all} ${cunits}\"\n"
            "print \"c23 = ${C23all} ${cunits}\"\n"
            "print \"c44 = ${C44all} ${cunits}\"\n"
            "print \"c55 = ${C55all} ${cunits}\"\n"
            "print \"c66 = ${C66all} ${cunits}\"\n"
            "print \"c14 = ${C14all} ${cunits}\"\n"
            "print \"c15 = ${C15all} ${cunits}\"\n"
            "print \"c16 = ${C16all} ${cunits}\"\n"
            "print \"c24 = ${C24all} ${cunits}\"\n"
            "print \"c25 = ${C25all} ${cunits}\"\n"
            "print \"c26 = ${C26all} ${cunits}\"\n"
            "print \"c34 = ${C34all} ${cunits}\"\n"
            "print \"c35 = ${C35all} ${cunits}\"\n"
            "print \"c36 = ${C36all} ${cunits}\"\n"
            "print \"c45 = ${C45all} ${cunits}\"\n"
            "print \"c46 = ${C46all} ${cunits}\"\n"
            "print \"c56 = ${C56all} ${cunits}\"\n"
            "print \"lammps_sim_done\"\n")
        return str_out

    def lammps_init_mod_to_string(self):
        atom_style = None
        if self.potential.is_charge:
            atom_style = 'charge'
        else:
            atom_style = 'atomic'
        str_out = (\
                "# ---- init.mod file\n"
                "variable up equal 1.0e-6\n"
                "units metal\n"
                "dimension 3\n"
                "boundary p p p\n"
                "atom_style {atom_style}\n"
                "atom_modify map array\n"
                "variable cfac equal 1.0e-4\n"
                "variable cunits string GPa\n"
                "# ---- define minimization parameters\n"
                "variable etol equal 0.0\n"
                "variable ftol equal 1.0e-10\n"
                "variable maxiter equal 100\n"
                "variable maxeval equal 1000\n"
                "variable dmax equal 1.0e-2\n"
                "# --- read data structure\n"
                "read_data structure.dat\n"
            ).format(\
                atom_style=atom_style
            )
        return str_out

    def lammps_displace_mod_to_string(self):
        str_out = (\
            "# NOTE: This script should not need to be\n"
            "# modified. See in.elastic for more info.\n"
            "# Find which reference length to use\n"
            "\n"
            "if \"${dir} == 1\" then &\n"
            "   \"variable len0 equal ${lx0}\"\n" 
            "if \"${dir} == 2\" then &\n"
            "   \"variable len0 equal ${ly0}\"\n" 
            "if \"${dir} == 3\" then &\n"
            "   \"variable len0 equal ${lz0}\"\n" 
            "if \"${dir} == 4\" then &\n"
            "   \"variable len0 equal ${lz0}\"\n" 
            "if \"${dir} == 5\" then &\n"
            "   \"variable len0 equal ${lz0}\"\n" 
            "if \"${dir} == 6\" then &\n"
            "   \"variable len0 equal ${ly0}\"\n" 
            "\n"
            "# Reset box and simulation parameters\n"
            "\n"
            "clear\n"
            "read_restart restart.equil remap\n"
            "include potential.mod\n"
            "\n"
            "# Negative deformation\n"
            "\n"
            "variable delta equal -${up}*${len0}\n"
            "if \"${dir} == 1\" then &\n"
            "   \"change_box all x delta 0 ${delta} remap units box\"\n"
            "if \"${dir} == 2\" then &\n"
            "   \"change_box all y delta 0 ${delta} remap units box\"\n"
            "if \"${dir} == 3\" then &\n"
            "   \"change_box all z delta 0 ${delta} remap units box\"\n"
            "if \"${dir} == 4\" then &\n"
            "   \"change_box all yz delta ${delta} remap units box\"\n"
            "if \"${dir} == 5\" then &\n"
            "   \"change_box all xz delta ${delta} remap units box\"\n"
            "if \"${dir} == 6\" then &\n"
            "   \"change_box all xy delta ${delta} remap units box\"\n"
            "\n"
            "# Relax atoms positions\n"
            "\n"
            "minimize ${etol} ${ftol} ${maxiter} ${maxeval}\n"
            "\n"
            "# Obtain new stress tensor\n"
            "\n"
            "variable tmp equal pxx\n"
            "variable pxx1 equal ${tmp}\n"
            "variable tmp equal pyy\n"
            "variable pyy1 equal ${tmp}\n"
            "variable tmp equal pzz\n"
            "variable pzz1 equal ${tmp}\n"
            "variable tmp equal pxy\n"
            "variable pxy1 equal ${tmp}\n"
            "variable tmp equal pxz\n"
            "variable pxz1 equal ${tmp}\n"
            "variable tmp equal pyz\n"
            "variable pyz1 equal ${tmp}\n"
            "\n"
            "# Compute elastic constant from pressure tensor\n"
            "\n"
            "variable C1neg equal ${d1}\n"
            "variable C2neg equal ${d2}\n"
            "variable C3neg equal ${d3}\n"
            "variable C4neg equal ${d4}\n"
            "variable C5neg equal ${d5}\n"
            "variable C6neg equal ${d6}\n"
            "\n"
            "# Reset box and simulation parameters\n"
            "\n"
            "clear\n"
            "read_restart restart.equil remap\n"
            "include potential.mod\n"
            "\n"
            "# Positive deformation\n"
            "\n"
            "variable delta equal ${up}*${len0}\n"
            "if \"${dir} == 1\" then &\n"
            "   \"change_box all x delta 0 ${delta} remap units box\"\n"
            "if \"${dir} == 2\" then &\n"
            "   \"change_box all y delta 0 ${delta} remap units box\"\n"
            "if \"${dir} == 3\" then &\n"
            "   \"change_box all z delta 0 ${delta} remap units box\"\n"
            "if \"${dir} == 4\" then &\n"
            "   \"change_box all yz delta ${delta} remap units box\"\n"
            "if \"${dir} == 5\" then &\n"
            "   \"change_box all xz delta ${delta} remap units box\"\n"
            "if \"${dir} == 6\" then &\n"
            "   \"change_box all xy delta ${delta} remap units box\"\n"
            "\n"
            "# Relax atoms positions\n"
            "\n"
            "minimize ${etol} ${ftol} ${maxiter} ${maxeval}\n"
            "\n"
            "# Obtain new stress tensor\n"
            "\n"
            "variable tmp equal pe\n"
            "variable e1 equal ${tmp}\n"
            "variable tmp equal press\n"
            "variable p1 equal ${tmp}\n"
            "variable tmp equal pxx\n"
            "variable pxx1 equal ${tmp}\n"
            "variable tmp equal pyy\n"
            "variable pyy1 equal ${tmp}\n"
            "variable tmp equal pzz\n"
            "variable pzz1 equal ${tmp}\n"
            "variable tmp equal pxy\n"
            "variable pxy1 equal ${tmp}\n"
            "variable tmp equal pxz\n"
            "variable pxz1 equal ${tmp}\n"
            "variable tmp equal pyz\n"
            "variable pyz1 equal ${tmp}\n"
            "\n"
            "# Compute elastic constant from pressure tensor\n"
            "\n"
            "variable C1pos equal ${d1}\n"
            "variable C2pos equal ${d2}\n"
            "variable C3pos equal ${d3}\n"
            "variable C4pos equal ${d4}\n"
            "variable C5pos equal ${d5}\n"
            "variable C6pos equal ${d6}\n"
            "\n"
            "# Combine positive and negative\n"
            "\n"
            "variable C1${dir} equal 0.5*(${C1neg}+${C1pos})\n"
            "variable C2${dir} equal 0.5*(${C2neg}+${C2pos})\n"
            "variable C3${dir} equal 0.5*(${C3neg}+${C3pos})\n"
            "variable C4${dir} equal 0.5*(${C4neg}+${C4pos})\n"
            "variable C5${dir} equal 0.5*(${C5neg}+${C5pos})\n"
            "variable C6${dir} equal 0.5*(${C6neg}+${C6pos})\n"
            "\n"
            "# Delete dir to make sure it is not reused\n"
            "\n"
            "variable dir delete\n")
        return str_out

    def postprocess(self):
        #LammpsSimulation.postprocess(self)
        self._get_elastic_tensor()

    def _get_elastic_tensor(self):
        filename = os.path.join(self.task_directory,'lammps.out')
        with open(filename,'r') as f:
            lines = f.readlines()

        lammps_results_names = ['c11','c12','c13','c14','c15','c16',
                                'c22','c23','c24','c25','c26',
                                'c33','c34','c35','c36',
                                'c44','c45','c46',
                                'c55','c56',
                                'c66']

        self.results = {}
        for i,line in enumerate(lines):
            for name in lammps_results_names:
                if '{} = '.format(name) in line:
                    try:
                        self.results[name] = \
                                float(line.split('=')[1].split()[0].strip())
                    except:
                        print('name:{}'.format(name))
                        print('line:{}'.format(line.strip()))
                        raise

        self.results['elastic_tensor'] = np.zeros((6,6))
        for i in range(6):
            for j in range(6):
                if i <= j:
                    self.results['elastic_tensor'][i,j] = \
                            self.results['c{}{}'.format(i+1,j+1)]
                else:
                    self.results['elastic_tensor'][i,j] = \
                            self.results['c{}{}'.format(j+1,i+1)]


                
