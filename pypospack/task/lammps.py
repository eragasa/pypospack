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
import pypospack.potential as potential
import pypospack.io.vasp as vasp
import pypospack.io.lammps as lammps
from pypospack.task import Task

atom_style_list = ['charge','atomic']

potential_map = {\
        'buckingham':['pypospack.potential','Buckingham'],
        'eam':['pypospack.potential','EmbeddedAtomModel'],
        'tersoff':['pypospack.potential','Tersoff']}

class LammpsSimulation(Task):
    """  Abstract data class for LAMMPS simulation objects 

    This is an abstract data class which defines the attributes and methods 
    necessary to interact with the Workflow manager.  A default implementation
    has been created for this class

    Args:
        task_name(str): unique id for the task name being defined
        task_directory(str): the directory where this task will create
            input and output files for LAMMPS.

    Attributes:
        config(:obj:'list' of :obj:'str'): a list of attributes required to 
            configure this LAMMPS task.
        config_map(dict):
        potential_map(dict):
    """
    def __init__(self,task_name,task_directory):
        self.status = None
        self.task_directory = ''
        self.structure_filename = ''
        self.config_dict = {\
                'potential_type':self.set_potential_type,
                'chemical_symbols':self.set_chemical_symbols,
                'structure_name':self.set_structure_name,
                'structure_filename':self.set_structure_filename}
        self.ready_dict = {\
                'param_dict':self.set_param_dict
                }
        self.run_dict = {}
        self.post_dict = {}
        self.done_dict = {}
        Task.__init__(self,task_name,task_directory)
        self.potential_map = potential_map
        if os.path.exists(self.task_directory):
            pass
        else:
            os.mkdir(self.task_directory)
        self.status='INIT'

    def run(self, is_mpi=False):
        self.write_lammps_input_file()
        if self.write_potential_files is None:
            self.write_potential_files(potential_files=self.potential_files)
        else:
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

    def lammps_input_file_to_string(self):
        str_out = "".join([\
                self._lammps_input_initialization_section(),
                self._lammps_input_create_atoms(),
                self._lammps_input_define_potential(),
                self._lammps_input_run_minimization(),
                self._lammps_input_out_section()])
        return(str_out)

    def write_lammps_input_file(self,filename='lammps.in'):
        str_out = self.lammps_input_file_to_string()
        filename = os.path.join(self.task_directory,filename)
        with open(filename,'w') as f:
            try:
                f.write(str_out)
            except:
                print('str_out:{}'.format(str_out))
                raise

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
    
        # read in poscar file
        vasp_structure = vasp.Poscar()
        try:
            vasp_structure.read(self.structure_filename)
        except filenotfounderror as e:
            print('error:',str(e))
            print('\tcurrent_working_directory',os.getcwd())
            print('\tself.structure_filename:',self.structure_filename)
            raise
        # instatiate using lammpsstructure file
        lammps_structure = lammps.LammpsStructure(\
                obj=vasp_structure)
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

        self.done_dict = {}
        for i,line in enumerate(lines):
            for name in lammps_result_names:
                if '{} = '.format(name) in line:
                    try:
                        self.done_dict[name] = \
                                float(line.split('=')[1].strip())
                    except:
                        print('name:{}'.format(name))
                        print('line:{}'.format(line.strip()))
                        raise

        self.done_dict['ecoh'] = self.done_dict['tot_energy'] / self.done_dict['num_atoms']

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
    def set_structure_name(self,name):
        self.structure_name =name

    def set_structure_filename(self,name):
        self.structure_filename = name

    def set_potential_type(self,name):
        self.potential_type = name

    def set_chemical_symbols(self,symbols):
        self.symbols = symbols

    def set_param_dict(self,param_dict):
        self.param_dict = copy.deepcopy(param_dict)

    def configure(self):
        self.configure_potential(self.potential_type)

    def configure_potential(self,pot_name):
        module_name = self.potential_map[pot_name][0]
        class_name = self.potential_map[pot_name][1]

        try:
            module = importlib.import_module(module_name)
            klass = getattr(module,class_name)
            self.potential = klass(self.symbols)
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

class LammpsMinimizeStructure(LammpsSimulation):
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

    def run(self,is_mpi=False):
        LammpsSimulation.run(self,is_mpi)

    def lammps_input_file_to_string(self):
        return LammpsSimulation.lammps_input_file_to_string(self)

    def write_lammps_input_file(self,filename='lammps.in'):
        LammpsSimulation.write_lammps_input_file(self,filename)

    def write_potential_files(self,param_dict=None,filename='potential.mod'):
        LammpsSimulation.write_potential_files(\
                self,param_dict,filename)
            
    def write_structure_file(self,filename='structure.dat'):
        LammpsSimulation.write_structure_file(self,filename)

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
        for k,v in ready_dict.items():
            self.ready_dict[k](v)
        self.status = 'READY'

    def req_run(self):
        """ this should method should be overridden, not inherited """
        return list(self.run_dict.keys())

    def send_run(self,run_dict):
        """ this should method should be overridden, not inherited """
        for k,v in run_dict.items():
            self.run_dict[k](v)
        self.run()
        self.status = "RUN"

    def req_post(self):
        """ this should method should be overridden, not inherited """
        return list(self.post_dict.keys())

    def send_post(self,post_dict):
        """ this should method should be overridden, not inherited """
        self.postprocess()
        self.status = "POST"

    def req_done(self):
        """ this should method should be overridden, not inherited """
        self.status = 'DONE'
        return list(self.done_dict)
        

    def send_done(self,done_dict):
        """ this should method should be overridden, not inherited """
        self.status = "DONE"

    def postprocess(self):
        LammpsSimulation.postprocess(self)
        if False:
            filename = os.path.join(self.task_directory,'lammps.out')
            with open(filename,'r') as f:
                lines = f.readlines()

            lammps_result_names = ['tot_energy','num_atoms',
                            'xx','yy','zz','xy','xz','yz',
                            'tot_press','pxx','pyy','pzz','pxy','pxz','pyz']

            self.done_dict = {}
            for i,line in enumerate(lines):
                for name in lammps_result_names:
                    if '{} = '.format(name) in line:
                        try:
                            self.done_dict[name] = \
                                    float(line.split('=')[1].strip())
                        except:
                            print('name:{}'.format(name))
                            print('line:{}'.format(line.strip()))
                            raise

            self.done_dict['ecoh'] = self.done_dict['tot_energy'] / self.done_dict['num_atoms']

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

