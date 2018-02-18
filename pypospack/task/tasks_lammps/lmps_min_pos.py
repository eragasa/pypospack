import os,copy
from collections import OrderedDict
from pypospack.io.vasp import Poscar
from pypospack.task.lammps import LammpsSimulation

class LammpsPositionMinimization(LammpsSimulation):
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
    def __init__(self,
            task_name,
            task_directory,
            structure_filename,
            restart=False,
            fullauto=False):
        _task_type = 'lmps_min_pos'

        self.bulk_structure_name = None
        self.bulk_structure_filename = None
        self.bulk_structure_lattice = None
        LammpsSimulation.__init__(self,
                task_name=task_name,
                task_directory=task_directory,
                task_type=_task_type,
                structure_filename=structure_filename,
                restart=restart,
                fullauto=fullauto)
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

    def on_init(self,configuration=None,results=None):
        LammpsSimulation.on_init(self,configuration=configuration)
        
        if 'bulk_structure' in configuration:
            self.bulk_structure_name = configuration['bulk_structure']
            self.bulk_structure_filename = configuration['bulk_structure_filename']
            self.bulk_structure_lattice = OrderedDict()

            _lattice_parameter_variables = [
                    'lmps_min_all.a11',
                    'lmps_min_all.a12',
                    'lmps_min_all.a13',
                    'lmps_min_all.a21',
                    'lmps_min_all.a22',
                    'lmps_min_all.a23',
                    'lmps_min_all.a31',
                    'lmps_min_all.a32',
                    'lmps_min_all.a33']

            self.bulk_lattice_components = OrderedDict()
            for k in _lattice_parameter_variables:
                _k = '{}.{}'.format(self.bulk_structure_name,k)
                self.bulk_lattice_components[_k] = None
    
    def on_config(self,configuration,results=None):
        if self.bulk_structure_name is not None:
            for k,v in results.items():
                if k in self.bulk_lattice_components:
                    self.bulk_lattice_components[k] = v
        LammpsSimulation.on_config(self,configuration=None,results=None)
    
    def get_conditions_ready(self):
        LammpsSimulation.get_conditions_ready(self)
        
        if self.bulk_structure_name is not None:
            _is_components_exist = []
            for k,v in self.bulk_lattice_components.items():
                _is_components_exist.append(v is not None)
            _all_components_exist = all(_is_components_exist)
            self.conditions_READY['bulk_lattice_components'] =\
                    _all_components_exist
    
    def on_ready(self,configuration=None,results=None):
        if self.bulk_structure_name is not None:
            self.__modify_structure(results=results)
        LammpsSimulation.on_ready(self,configuration=configuration)

    def __modify_structure(self,results):
        assert isinstance(results,dict)
        #_ideal_structure_name = self.ideal_structure_name
        a11_n = '{}.lmps_min_all.a11'.format(self.bulk_structure_name)
        a12_n = '{}.lmps_min_all.a12'.format(self.bulk_structure_name)
        a13_n = '{}.lmps_min_all.a13'.format(self.bulk_structure_name)
        a21_n = '{}.lmps_min_all.a21'.format(self.bulk_structure_name)
        a22_n = '{}.lmps_min_all.a22'.format(self.bulk_structure_name)
        a23_n = '{}.lmps_min_all.a23'.format(self.bulk_structure_name)
        a31_n = '{}.lmps_min_all.a31'.format(self.bulk_structure_name)
        a32_n = '{}.lmps_min_all.a32'.format(self.bulk_structure_name)
        a33_n = '{}.lmps_min_all.a33'.format(self.bulk_structure_name)
        a11 = self.bulk_lattice_components[a11_n]
        a12 = self.bulk_lattice_components[a12_n]
        a13 = self.bulk_lattice_components[a13_n]
        a0 = (a11*a11+a12*a12+a13*a13)**0.5
        self.structure.a0 = a0

        
    def on_post(self,configuration=None):
        self.__get_results_from_lammps_outputfile()
        LammpsSimulation.on_post(self,configuration=configuration)
    
    def __get_results_from_lammps_outputfile(self):
        _filename = os.path.join(
                self.task_directory,
                'lammps.out')
        with open(_filename,'r') as f:
            lines = f.readlines()
        
        _variables = [
                'tot_energy',
                'num_atoms',
                'xx','yy','zz','xy','xz','yz',
                'tot_press',
                'pxx', 'pyy', 'pzz', 'pxy', 'pxz', 'pyz',
                ]
        _results = OrderedDict()
        
        for i,line in enumerate(lines):
            for name in _variables:
                if line.startswith('{} = '.format(name)):
                    _results[name] = float(line.split('=')[1].strip())

                if line.startswith('ERROR:'):
                    print('name:{}'.format(name))
                    print('line:{}'.format(line.strip))
                    raise NotImplementedError
      
        _task_name = self.task_name
        self.results = OrderedDict()
        self.results['{}.{}'.format(_task_name,'toten')] = _results['tot_energy']
        self.results['{}.{}'.format(_task_name,'natoms')] = _results['num_atoms']
        # this only works for orthogonal cells
        self.results['{}.{}'.format(_task_name,'a11')] = _results['xx']
        self.results['{}.{}'.format(_task_name,'a12')] = 0
        self.results['{}.{}'.format(_task_name,'a13')] = 0
        self.results['{}.{}'.format(_task_name,'a21')] = 0
        self.results['{}.{}'.format(_task_name,'a22')] = _results['yy']
        self.results['{}.{}'.format(_task_name,'a23')] = 0
        self.results['{}.{}'.format(_task_name,'a31')] = 0
        self.results['{}.{}'.format(_task_name,'a32')] = 0
        self.results['{}.{}'.format(_task_name,'a33')] = _results['zz']
        self.results['{}.{}'.format(_task_name,'totpress')] = _results['tot_press']
        self.results['{}.{}'.format(_task_name,'p11')] = _results['pxx']
        self.results['{}.{}'.format(_task_name,'p12')] = _results['pxy']
        self.results['{}.{}'.format(_task_name,'p13')] = _results['pxz']
        self.results['{}.{}'.format(_task_name,'p21')] = _results['pxy']
        self.results['{}.{}'.format(_task_name,'p22')] = _results['pyy']
        self.results['{}.{}'.format(_task_name,'p23')] = _results['pyz'] #pyz=pzy
        self.results['{}.{}'.format(_task_name,'p31')] = _results['pxz'] #pxz=pzx
        self.results['{}.{}'.format(_task_name,'p32')] = _results['pyz']
        self.results['{}.{}'.format(_task_name,'p33')] = _results['pzz']

    def _lammps_input_run_minimization(self):
        str_out = (
            '# ---- define settings\n'
            'compute eng all pe/atom\n'
            'compute eatoms all reduce sum c_eng\n'
            '\n'
            '# ---- run minimization\n'
            'reset_timestep 0\n'
            'thermo 1\n'
            'thermo_style custom step pe lx ly lz xy xz yz press pxx pyy pzz pxy pxz pyz c_eatoms\n'
            'min_style cg\n'
            'minimize 1e-20 1e-20 1000 100000\n'
            '\n'
            )
        return str_out

