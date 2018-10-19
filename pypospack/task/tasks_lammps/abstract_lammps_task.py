""" Implementation of AbstractLammpsSimulation

"""
import os,copy,importlib,subprocess
from collections import OrderedDict

import pypospack.io.vasp as vasp
import pypospack.io.lammps as lammps

from pypospack.task import Task

from pypospack.io.eamtools import EamSetflFile
from pypospack.exceptions import LammpsSimulationError
from pypospack.potential import Potential,EamPotential,PotentialObjectMap
from pypospack.potential import StillingerWeberPotential

class AbstractLammpsSimulation(Task):
    """ Calculates cohesive energy

    This is an abstract data class which defines the attributes and methods 
    necessary to interact with the Workflow manager.  A default implementation
    has been created for this class.  This class calculates the energy of
    a simulation cell.  No structural or positions are calculated.

    Args:
        task_name(str): unique id for the task name being defined
        task_directory(str): the directory where this task will create
            input and output files for LAMMPS.
        structure_filename(str):
        simulation_type(str):
        restart(bool)
        fullauto(bool)
    Attributes:
        task_name(str)
        task_directory(str)
        task_type(str)
        is_restart(bool)
        is_fullauto(bool)
        lammps_input_filename(str)
        lammps_output_filename(str)
        lammps_structure_filename(str)
        lammps_eam_filename(str)
        potential(pypospack.potential.Potential): the potential class
        structure_filename(str)
        structure(pypospack.io.vasp.Poscar): the structure
        config(:obj:'list' of :obj:'str'): a list of attributes required to 
            configure this LAMMPS task.
        config_map(dict):
        potential_map(dict):
        results(dict): results of the simulation
        lammps_bin(str): location of the serial lammps binary
        conditions_INIT(collections.OrderedDict)
        conditions_CONFIG(collections.OrderedDict)
        conditions_READY(collections.OrderedDict)
        conditions_RUNNING(collections.OrderedDict)
        conditions_POST(collections.OrderedDict)
        conditions_FINISHED(collections.OrderedDict)
        conditions_ERROR(collections.OrderedDict)
    """
    def __init__(self,
            task_name,
            task_directory,
            task_type='single_point',
            task_requires=None,
            structure_filename='POSCAR',
            restart=False,
            fullauto=False):

        self.is_restart = restart
        self.is_fullauto = fullauto

        self.potential = None
        self.structure = None
        self.structure_filename = structure_filename
        self.structure = vasp.Poscar()
        self.structure.read(self.structure_filename)

        self.lammps_input_filename = 'lammps.in'
        self.lammps_output_filename = 'lammps.out'
        self.lammps_structure_filename = 'lammps.structure'
        self.lammps_potentialmod_filename = 'potential.mod'
        self.lammps_setfl_filename = None
        self.lammps_bin = os.environ['LAMMPS_BIN']

        # flowcontrol filename
        self.results = None
        self.results_filename = 'pypospack.{}.out'.format(task_name)

        self.process = None
        self.results_processed = None
        # configuration
        self.configuration = OrderedDict()
        self.task_type = task_type
        
        self.conditions_INIT = None
        self.conditions_CONFIG = None
        self.conditions_READY = None
        self.conditions_RUNNING = None
        self.conditions_POST = None
        self.conditions_FINISHED = None
        self.conditions_ERROR = None
        Task.__init__(
                self,
                task_name=task_name,
                task_directory=task_directory,
                restart=restart)
        self.task_requires = copy.deepcopy(task_requires)

    
    @property
    def parameters(self):
        if self.potential is None:
            return None
        else:
            return self.potential.parameters

    @parameters.setter
    def parameters(self,parameters):
        self.potential.parameters = parameters
   
    @property
    def symbols(self):
        return self.potential.symbols

    def set_potential_parameters(self,parameters=None):
        """
        Sets the parameters of the potential

        Args:
            parameters(OrderedDict)
        """
        if self.potential is None:
            return
        
        if parameters is not None:
            _parameters = copy.deepcopy(parameters)
            self.parameters = _parameters

        if parameters in self.configuration:
            _parameters = copy.deepcopy(self.configuration['parameters'])
            self.parameters = _parameters
    
    def on_init(self,configuration=None):
        if configuration is not None:
            self.configuration = copy.deepcopy(configuration)
       
        try:
            self.configure_potential(potential=self.configuration['potential'])
        except KeyError:
            print("self.configuration:",type(self.configuration))
            for k,v in self.configuration.items():
                print(k,' = ',v)
            raise

        # writing eam potential files
        if type(self.potential) is EamPotential:
            if self.lammps_setfl_filename is None:
                self.lammps_setfl_filename = "{}.eam.alloy".format(
                        "".join(self.potential.symbols))

            # if setfl_filename_src is set, then we just copy the 
            # EAM potential file.

            
            if all([self.potential.obj_pair is None,
                    self.potential.obj_density is None,
                    self.potential.obj_embedding is None,
                    self.potential.setfl_filename_src is not None]):
                _eam_setfl_filename_src = self.potential.setfl_filename_src
                _eam_setfl_filename_dst = os.path.join(
                        self.task_directory,
                        self.lammps_setfl_filename)
                shutil.copyfile(
                        src=_eam_setfl_filename_src,
                        dst=_eam_setfl_filename_dst)
            elif all([self.potential.obj_pair is not None,
                      self.potential.obj_density is not None,
                      self.potential.obj_embedding is not None]):
                pass
            else:
                msg_err = (
                    "EamPotential must be either be parameterized by setting "
                    "pair,density,and embedding formalisms through the "
                    "constructor or a setfl filename must be provided through "
                    "the filename argument\n"
                    "obj_pair:{obj_pair}\n"
                    "obj_density:{obj_pensity}\n"
                    "obj_embedding:{obj_embedding}\n"
                    "setfl_filename:{setfl_filename\n"
                    ).format(
                            obj_pair=str(type(self.potential.obj_pair)),
                            obj_density=str(type(self.potential.obj_density)),
                            obj_embedding=str(type(self.potential.obj_embedding)),
                            setfl_filename=str(self.potential.setfl_filename))
                raise ValueError(msg_err)
        
        self.set_potential_parameters()
        if self.configuration['parameters'] is not None:
            self.parameters = self.configuration['parameters']

        if self.structure is None:
            if self.structure_filename is not None:
                self.read_structure_file()

        self.update_status()
        if self.is_fullauto:
            self.on_update_status()


    def on_config(self,configuration=None,results=None):
        if configuration is not None:
            self.configuration = copy.deepcopy(configuration)

        self.configure_potential()
        if 'parameters' in self.configuration:
            if isinstance(self.potential,Potential):
                _parameters = self.configuration['parameters']
                self.potential.parameters = _parameters
        
        # writing eam potential files
        if type(self.potential) is EamPotential:
            if self.lammps_setfl_filename is None:
                self.lammps_setfl_filename = "{}.eam.alloy".format(
                        "".join(self.potential.symbols))

            # if setfl_filename_src is set, then we just copy the 
            # EAM potential file.
            if all([self.potential.obj_pair is None,
                    self.potential.obj_density is None,
                    self.potential.obj_embedding is None,
                    self.potential.setfl_filename_src is not None]):
                _eam_setfl_filename_src = self.potential.setfl_filename_src
                _eam_setfl_filename_dst = os.path.join(
                        self.task_directory,
                        self.lammps_setfl_filename)
                shutil.copyfile(
                        src=_eam_setfl_filename_src,
                        dst=_eam_setfl_filename_dst)
            elif all([self.potential.obj_pair is not None,
                      self.potential.obj_density is not None,
                      self.potential.obj_embedding is not None]):
                pass
            else:
                msg_err = (
                    "EamPotential must be either be parameterized by setting "
                    "pair,density,and embedding formalisms through the "
                    "constructor or a setfl filename must be provided through "
                    "the filename argument\n"
                    "obj_pair:{obj_pair}\n"
                    "obj_density:{obj_pensity}\n"
                    "obj_embedding:{obj_embedding}\n"
                    "setfl_filename:{setfl_filename\n"
                    ).format(
                            obj_pair=str(type(self.potential.obj_pair)),
                            obj_density=str(type(self.potential.obj_density)),
                            obj_embedding=str(type(self.potential.obj_embedding)),
                            setfl_filename=str(self.potential.setfl_filename))
                raise ValueError(msg_err)

        self.update_status()
        if self.is_fullauto:
            self.on_update_status()

    def on_ready(self,configuration=None,results=None):
        if configuration is not None:
            self.configuration = copy.deepcopy(configuration)

        self.write_lammps_input_file()
        self.write_potential_file()
        self.write_structure_file()
        if isinstance(self.potential,EamPotential):
            if self.potential.setfl_filename_src is None:
                if self.lammps_setfl_filename is None:
                    self.lammps_setfl_filename = '{}.eam.alloy'.format(
                            "".join(self.potential.symbols))
                _eam_setfl_filename_dst = os.path.join(
                        self.task_directory,
                        self.lammps_setfl_filename)
                self.write_eam_setfl_file(
                        filename=_eam_setfl_filename_dst)
        self.run()

        self.update_status()
        if self.is_fullauto:
            self.on_update_status()

    def on_running(self,configuration=None):
        self.update_status()
        if self.is_fullauto:
            self.on_update_status()
   
    def on_post(self,configuration=None):
        self.results_processed = True
        self.update_status()
        if self.is_fullauto:
            self.on_update_status()

    def on_finished(self,configuration=None):
        # doing nothing here
        self.update_status()
        if self.is_fullauto:
            self.on_update_status()

    def on_error(self,configuration=None):
        raise ValueError()

    def get_conditions_init(self):
        self.conditions_INIT = OrderedDict()
        self.conditions_INIT['task_directory_created']\
                = os.path.isdir(self.task_directory)
        return self.conditions_INIT

    def get_conditions_config(self):
        self.conditions_CONFIG = OrderedDict()
        self.conditions_CONFIG['potential_initialized']\
                = isinstance(self.potential, Potential)

        if type(self.potential) is EamPotential \
                and self.potential.obj_pair is None \
                and self.potential.obj_density is None \
                and self.potential.obj_embedding is None \
                and self.potential.setfl_filename_src is not None:
            self.conditions_CONFIG['parameters_processed'] = True
            self.conditions_CONFIG['potential_is_eam_setfl_file'] = True
        elif self.potential is None:
            self.conditions_CONFIG['parameters_processed'] = False
        else:
            self.conditions_CONFIG['parameters_processed']\
                    = all([v is not None 
                        for k,v in self.potential.parameters.items()])

    def get_conditions_ready(self):
        self.conditions_READY = OrderedDict()
    
    def get_conditions_running(self):
        self.conditions_RUNNING = OrderedDict()
        self.conditions_RUNNING['process_initialized'] = self.process is not None
   
    def get_conditions_post(self):
        self.conditions_POST = OrderedDict()
        
        if self.process is None:
            _process_finished = False
        else:
            _poll_result = self.process.poll()
            if _poll_result is None:
                _process_finished = False
            elif _poll_result == 0:
                _process_finished = True
            elif _poll_results == 1:
                if self.conditions_ERROR is None:
                    self.conditions_ERROR=OrderedDict()

                m = "Lammps excited with status {}.  If running an EAM "
                m += "potential this is most likely caused by an out-of-index "
                m += "exception because the electron density is too high when "
                m += "evaluating the embedding function.  The code for modifying "
                m += "max_rho for the embedding function is in pypospack.potential.EamPotential"

                self.conditions_ERROR['lmps_bin_err') = err_msg
                raise LammpsSimulationError(err_msg)
            else:
                if self.conditions_ERROR is None:
                    self.conditions_ERROR= OrderedDict()
                err_msg = 'Lammps exited with status {}.'.format(_poll_result)
                self.conditions_ERROR['lmps_bin_err'] =  err_msg
                raise LammpsSimulationError(err_msg)
        
        self.conditions_POST['process_finished'] = _process_finished


    def get_conditions_finished(self):
        self.conditions_FINISHED = OrderedDict()
        self.conditions_FINISHED['is_processed'] = self.results_processed  
    
    def get_conditions_error(self):
        if self.conditions_ERROR is None:
            self.conditions_ERROR = OrderedDict()

    def restart(self):
        raise NotImplementedError

    def run(self):

        _lammps_bin = self.lammps_bin
        cmd_str = '{} -i lammps.in > lammps.out'.format(_lammps_bin)

        # change context directory
        
        _cwd = os.getcwd()
        os.chdir(self.task_directory)
        
        # https://stackoverflow.com/questions/4789837/how-to-terminate-a-python-subprocess-launched-with-shell-true/4791612#4791612
        self.process = subprocess.Popen(
                cmd_str,
                shell=True,
                cwd=self.task_directory,
                preexec_fn=os.getpgrp)
                #preexec_fn=os.setsid)
        self.process_info = OrderedDict()
        os.chdir(_cwd)

    def post(self):
        lammps_result_names = ['tot_energy','num_atoms',
                        'xx','yy','zz','xy','xz','yz',
                        'tot_press','pxx','pyy','pzz','pxy','pxz','pyz']

        self.results = OrderedDict()
        self.get_variables_from_lammps_output(
                variables = lammps_result_names)
        
        try:
            # calculate cohesive energy
            total_energy = self.results['tot_energy']
            n_atoms = self.results['num_atoms']
            self.results['ecoh'] = total_energy/n_atoms
        except KeyError as e:
            print(e)

        self.status = 'DONE'

    def get_variables_from_lammps_output(self,variables):
        filename = os.path.join(self.task_directory,'lammps.out')
        with open(filename,'r') as f:
            lines = f.readlines()

        self.results = {}
        for i,line in enumerate(lines):
            for name in variables:
                if line.startswith('{} = '.format(name)):
                    try:
                        self.results[name] = \
                                float(line.split('=')[1].strip())
                    except ValueError as e:
                        if line.split('=')[1].strip().endswith('GPa'):
                            self.results[name] = \
                                float(line.split('=')[1].strip().split(' ')[0])
                        else:
                            raise
                    except:
                        print('name:{}'.format(name))
                        print('line:{}'.format(line.strip()))
                        raise

    def write_eam_setfl_file(self,filename):
        _setfl_dst_filename = filename
        _Nr = self.configuration['potential']['N_r']
        _Nrho = self.configuration['potential']['N_rho']
        _rmax = self.configuration['potential']['r_max']
        _rhomax = self.configuration['potential']['rho_max']
        _rcut = self.configuration['potential']['r_cut']
        _parameters = self.configuration['parameters']

        _a0 = self.configuration['potential']['a0']
        _latt_type = self.configuration['potential']['lattice_type']

        _rcut = self.potential.determine_r_max(a0=_a0,latt_type='fcc')
        _rhomax = self.potential.determine_rho_max(a0=_a0,latt_type='fcc')
        self.configuration['potential']['rho_max'] = _rhomax
        self.configuration['potential']['rcut'] = _rcut

        self.potential.write_setfl_file(
                filename=_setfl_dst_filename,
                symbols=self.potential.symbols,
                Nr=_Nr,
                rmax=_rmax,
                rcut=_rcut,
                Nrho=_Nrho,
                rhomax=_rhomax,
                parameters=_parameters)

    def write_potential_file(self):
        if self.potential is None:
            return
        
        _setfl_dst_filename = None

        # <-------- FOR EAM POTENTIALS
        if isinstance(self.potential,EamPotential):
            _symbols = "".join(self.potential.symbols)
            _filename = "{}.eam.alloy".format(_symbols)
            _setfl_dst_filename = os.path.join(self.task_directory,_filename)
            _str_out = self.potential.lammps_potential_section_to_string(
                setfl_dst_filename=_setfl_dst_filename)
        
        # <-------- FOR STILLINGER WEBER POTENTIALS
        elif isinstance(self.potential,StillingerWeberPotential):
            # define the filename --- SiO.parameters, Si.parameters
            _symbols_str = "".join(self.potential.symbols)
            _p_fname = "{}.parameters".format(_symbols_str)

            # set the name of the output file
            self.potential.lmps_parameter_filename = _p_fname
            
            # get the string of potential.mod
            _str_out = self.potential.lammps_potential_section_to_string()

            # write out the potential parameter file
            _str_lmps_params = self.potential.lammps_parameter_file_to_string()
            
            _p_fname_dst = os.path.join(self.task_directory,_p_fname)
            with open(_p_fname_dst,'w') as f:
                f.write(_str_lmps_params)

        else:
            _str_out = self.potential.lammps_potential_section_to_string()

        _str_out += "\n"
        
        # <-------- EWALD CHARGE SUMMATION METHOD 
        if self.potential.is_charge:
            _str_out += "kspace_style pppm 1.0e-5\n"
            _str_out += "\n"
        
        # <-------- TREATMENT OF NEAREST NEIGHBORS
        _str_out += "neighbor 1.0 bin\n"
        _str_out += "neigh_modify every 1 delay 0 check yes\n"

        # <-------- WRITE POTENTIAL.MOD TO FILESYSTEM
        _lammps_potentialmod_filename = os.path.join(
                self.task_directory,
                self.lammps_potentialmod_filename)
        with open(_lammps_potentialmod_filename,'w') as f:
            f.write(_str_out)
    
    def write_lammps_input_file(self,filename='lammps.in'):
        """ writes LAMMPS input file 
        
        Args:
            filename (str): name of the input file for LAMMPS. Default is
                'lammps.in'.
        """
        str_out = self.lammps_input_file_to_string()
        filename = os.path.join(self.task_directory,filename)
        with open(filename,'w') as f:
            f.write(str_out)

    def lammps_input_file_to_string(self):
        """ string for the LAMMPS input file """

        str_out = "".join([\
                self._lammps_input_initialization_section(),
                self._lammps_input_create_atoms(),
                self._lammps_input_define_potential(),
                self._lammps_input_run_minimization(),
                self._lammps_input_out_section()])
        return(str_out)
            
    def write_structure_file(self,filename=None):
        if filename is not None:
            self.lammmps_structure_filename = filename
        
        _filename = os.path.join(
                self.task_directory,
                self.lammps_structure_filename)
        if self.potential.is_charge:
            _atom_style = 'charge'
        else:
            _atom_style = 'atomic'
        _symbol_list = self.potential.symbols
    
        # instatiate using lammpsstructure file
        self.lammps_structure = lammps.LammpsStructure(\
                obj=self.structure)
        
        self.lammps_structure.write(\
                filename=_filename,
                symbol_list=_symbol_list,
                atom_style=_atom_style)

    def modify_structure(self):pass

    def configure_potential(self,potential=None):
        """
            Args:
                potential(OrdereDict): potential is a dictionary which has the
                    necessary keywords to configure the object.

            Notes:
                For buckingham potential,
                    potential[potential_type] = 'buck'
                    potential['symbols'] = ['Mg','O']
                    potential['global_cutoff'] = [10.0]
                For eam potentials,
                    potential['potential_type'] = 'eam'
                    potential['pair_type'] = 'morse'
                    potential['density_type'] = 'eam_exp_dens'
                    potential['embedding_type'] = 'eam_universal_embedding'
        """
        assert isinstance(potential,OrderedDict) or potential is None

        _potential = None
        if isinstance(potential,OrderedDict):
            self.configuration['potential']=copy.deepcopy(potential)
            _potential = self.configuration['potential']
        elif potential is None:
            if isinstance(self.potential,Potential):
                return
            else:
                _potential = self.configuration['potential']

        _potential_type = _potential['potential_type']
        if _potential_type == 'eam':
            
            if 'setfl_filename' not in _potential:
                _potential['setfl_filename'] = None
                self.configuration['potential']['setfl_filename'] = None

            # configure the eam potential using an external file
            if _potential['setfl_filename'] is not None:
                _module_name,_class_name = PotentialObjectMap(
                        potential_type=_potential_type)
                _symbols = self.configuration['potential']['symbols']
                _setfl_filename = _potential['setfl_filename']
                try:
                    _module = importlib.import_module(_module_name)
                    _class = getattr(_module,_class_name)
                    self.potential = _class(
                            symbols=_symbols,
                            filename=_setfl_filename)
                except:
                    raise
            # configure the eam potential using parameters
            else:
                # get information from the configuration dictionary
                _module_name,_class_name = PotentialObjectMap(
                        potential_type=_potential_type)
                _symbols = self.configuration['potential']['symbols']
                _pair_type = _potential['pair_type']
                _embedding_type = _potential['embedding_type']
                _density_type = _potential['density_type']
                try:
                    _module = importlib.import_module(_module_name)
                    _class = getattr(_module,_class_name)
                    self.potential = _class(
                            symbols=_symbols,
                            func_pair=_pair_type,
                            func_embedding=_embedding_type,
                            func_density=_density_type)
                except:
                    raise
        # <-------- parameterized classes which don't have stupid tabularized
        #           methods for representing potentials, like EAM.
        else:
            _module_name,_class_name = PotentialObjectMap(
                    potential_type=_potential_type)
            _symbols = self.configuration['potential']['symbols']
            # all other potentials are parameterized
            try:
                _module = importlib.import_module(_module_name)
                _class = getattr(_module,_class_name)
                self.potential = _class(symbols=_symbols)
            except:
                raise

    # private functions for building lammps input files
    def _lammps_input_initialization_section(self):
        if self.potential.is_charge:
            _atom_style = 'charge'
        else:
            _atom_style = 'atomic'
        str_out = (
            '# ---- initialize simulations\n'
            'clear\n'
            'units metal\n'
            'dimension 3\n'
            'boundary p p p\n'
            'atom_style {atom_style}\n'
            'atom_modify map array\n'
            ).format(
                    atom_style=_atom_style)
        return str_out

    def _lammps_input_create_atoms(self):
        _structure_filename = os.path.join(
                self.task_directory,
                self.lammps_structure_filename)
        str_out = (
            '# ---- create atoms\n'
            'read_data {structure_filename}\n'
            ).format(
                    structure_filename=_structure_filename)

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
            'fix 1 all box/relax aniso 0.0 vmax 0.001\n'
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
            'variable a11 equal "xhi-xlo"\n'
            'variable a22 equal "yhi-ylo"\n'
            'variable a33 equal "zhi-zlo"\n'
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
            'print \"a11 = ${a11}\"\n'
            'print \"a22 = ${a22}\"\n'
            'print \"a33 = ${a33}\"\n'
            'print \"a12 = ${tilt_xy}\"\n'
            'print \"a13 = ${tilt_xz}\"\n'
            'print \"a23 = ${tilt_yz}\"\n'
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
