import os,copy
from collections import OrderedDict
from pypospack.task.lammps import LammpsSimulation

class LammpsStaticCalculations(LammpsSimulation):
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

        _task_type = 'lmps_min_none'
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

    def on_config(self,configuration,results=None):
        LammpsSimulation.on_config(self,configuration=None,results=None)
    
    def on_post(self,configuration=None):
        self.__get_results_from_lammps_outputfile()
        LammpsSimulation.on_post(self,configuration=configuration)
    
    def on_ready(self,configuration=None,results=None):
        LammpsSimulation.on_ready(
                self,
                configuration=configuration,
                results=results)
    
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
        self.results['{}.{}'.format(_task_name,'a11')] = _results['a11']
        self.results['{}.{}'.format(_task_name,'a12')] = _results['a12']
        self.results['{}.{}'.format(_task_name,'a13')] = _results['a13']
        self.results['{}.{}'.format(_task_name,'a21')] = 0
        self.results['{}.{}'.format(_task_name,'a22')] = _results['a22']
        self.results['{}.{}'.format(_task_name,'a23')] = _results['a23']
        self.results['{}.{}'.format(_task_name,'a31')] = 0
        self.results['{}.{}'.format(_task_name,'a32')] = 0
        self.results['{}.{}'.format(_task_name,'a33')] = _results['a33']
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
            '# ---- run minimization\n'            
            'reset_timestep 0\n'
            'thermo 10\n'
            'thermo_style custom step pe lx ly lz xy xz yz press pxx pyy pzz pxy pxz pyz c_eatoms\n'
            'run 0\n"
            )
        return str_out

