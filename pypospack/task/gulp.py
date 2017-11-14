# -*- coding: utf-8 -*-
""" Implementation of GulpSimulation and other implemented classes

This module implements simulation tasks to integrate with GULP.  A template
class is implemented as GulpSimulation, which should be subclassed for new
implementations requiring GULP simulations.
"""

import os, copy, importlib, subprocess
from collections import OrderedDict
from pathlib import Path
import pypospack.potential as potential
import pypospack.crystal as crystal
import pypospack.io.vasp as vasp
import pypospack.io.gulp as gulp
from pypospack.task import Task

potential_map = {\
        'buckingham':['pypospack.potential','Buckingham'],
        'eam':['pypospack.potential','EmbeddedAtomModel'],
        'tersoff':['pypospack.potential','Tersoff']}

class GulpSimulationError():
    pass

class GulpSimulation(Task):
    def __init__(self,
            task_name,
            task_directory,
            structure_filename='POSCAR',
            restart=False,
            fullauto=False):
      
        self.is_fullauto = fullauto
        self.running_filename = 'running'
        self.simulation_type = 'gulp'

        # gulp information
        self.gulp_input_filename = 'gulp.in'
        self.gulp_output_filename = 'gulp.out'
        self.gulp_bin = os.environ['GULP_BIN']
        
        # structure information
        self.structure_filename = structure_filename
        self.structure = vasp.Poscar()
        self.structure.read(self.structure_filename)

        # potential information
        self.potential = None
        self._parameters = None
        self.eam_setfl_filename = None

        # minimization flags
        self.is_optimize = True
        self.is_constant_pressure = True
        self.is_constant_volume = False
        
        # configuration
        self.configuration = OrderedDict()
        self.configuration['potential'] = None
        self.configuration['parameters'] = None
        
        # results
        self.results = OrderedDict()
        self.results_filename = '{}.results.dat'.format(task_name)
        # status
        # initialize the parent class
        Task.__init__(self,
                task_name=task_name,
                task_directory=task_directory,
                restart=restart)

    def on_init(self,configuration=None):
        if configuration is not None:
            self.configuration = copy.deepcopy(configuration)
        if self.configuration['potential'] is not None:
            if self.potential is None:
                self.configure_potential()
        if self.configuration['parameters'] is not None:
            self.parameters = self.configuration['parameters']
        if self.structure is None:
            if self.structure_filename is not None:
                self.read_structure_file()
        
        self.update_status()
        if self.is_fullauto:
            self.on_update_status()

    def on_config(self,configuration=None):
        if configuration is not None:
            self.configuration = copy.deepcopy(configuration)
        _gulp_input_filename = os.path.join(
                self.task_directory,
                self.gulp_input_filename)
        if not os.path.exists(_gulp_input_filename):
            self.write_gulp_input_file(
                    filename=_gulp_input_filename)

        self.update_status()
        if self.is_fullauto:
            self.on_update_status()
    def on_ready(self,configuration=None):

        if configuration is not None:
            self.configuration = copy.deepcopy(configuration)

        running_filename = os.path.join(
                self.task_directory,
                self.running_filename)
        Path(running_filename).touch()
        
        self.run()

        self.update_status()
        if self.is_fullauto:
            self.on_update_status()

    def on_running():
        pass

    def on_post():
         pass

    def on_error():
        pass

    # override default behavior
    def get_conditions_init(self):
        self.conditions_INIT = OrderedDict()
        self.conditions_INIT['is_task_directory_created'] \
                = os.path.isdir(self.task_directory)

    # override default behavior
    def get_conditions_config(self):
        self.conditions_CONFIG = OrderedDict()
        self.conditions_CONFIG['is_structure_configured'] \
                = isinstance(
                        self.structure,
                        crystal.SimulationCell)
        self.conditions_CONFIG['is_potential_configured'] \
                = isinstance(
                        self.potential,
                        potential.Potential)
        try: 
            self.conditions_CONFIG['has_all_parameters'] \
                    = all([v is not None for k,v in self.parameters.items()])
        except AttributeError as e:
            if str(e) == "'NoneType' object has no attribute 'items'":
                self.conditions_CONFIG['has_all_parameters'] = False

    def get_conditions_ready(self):
        self.conditions_READY = OrderedDict()
        self.conditions_READY['input_file_is_written'] \
                = os.path.isfile(os.path.join(
                    self.task_directory,
                    self.gulp_input_filename))
    
    def get_conditions_running(self):
        self.conditions_RUNNING = OrderedDict()
        self.conditions_RUNNING['is_running_file'] \
                = os.path.isfile(os.path.join(
                    self.task_directory,
                    self.running_filename))
    
    def get_conditions_post(self):
        self.conditions_POST = OrderedDict()
        self.conditions_POST['is_gulp_out_exists'] \
                = os.path.isfile(os.path.join(
                    self.task_directory,
                    self.running_filename))

    def get_conditions_finished(self):
        self.conditions_FINISHED = OrderedDict()
        self.conditions_FINISHED['is_results_exists'] \
                = os.path.isfile(self.results_filename)

    def get_conditions_error(self):
        self.conditions_ERROR = OrderedDict()


    def run(self):
        gulp_input_filename = os.path.join(
                self.task_directory,
                self.gulp_input_filename)
        if not os.path.isfile(gulp_input_filename):
            self.write_gulp_input_file(filename=gulp_input_filename)

        gulp_output_filename = os.path.join(
                self.task_directory,
                self.gulp_output_filename)
        if not os.path.isfile(gulp_output_filename):
            cmd_str = '{gulp_bin} < {gulp_in} > {gulp_out} 2>/dev/null'.format(
                    gulp_bin=self.gulp_bin,
                    gulp_in=self.gulp_input_filename,
                    gulp_out=self.gulp_output_filename)

            os.chdir(self.task_directory)
            try:
                subprocess.call(cmd_str,shell=True)
            finally:
                os.chdir(self.root_directory)

    def read_structure_file(self,filename=None):
        if filename is not None:
            self.structure_filename

        self.structure = vasp.Poscar()
        self.structure.read(self.structure_filename)
    
    def write_gulp_input_file(self,filename=None):
        if filename is not None:
            self.gulp_input_filename = 'gulp.in'

        str_gulp_in = "opti prop conp\n"
        str_gulp_in += self.get_gulpinputfile_structuresection()
        str_gulp_in += self.get_gulpinputfile_potentialsection()

        _filename = os.path.join(
                self.task_directory,
                self.gulp_input_filename)
        with open(_filename,'w') as f:
            f.write(str_gulp_in)

    def get_gulpinputfile_structuresection(self,filename=None):
        if filename is not None:
            self.structure_filename = filename

        # read in the poscar file
        self.structure = vasp.Poscar()
        self.structure.read(self.structure_filename)

        H = self.structure.H * self.structure.a0
        str_out = "vectors\n"
        str_out += "{:.10f} {:.10f} {:.10f}\n".format(H[0,0],H[0,1],H[0,2])
        str_out += "{:.10f} {:.10f} {:.10f}\n".format(H[1,0],H[1,1],H[1,2])
        str_out += "{:.10f} {:.10f} {:.10f}\n".format(H[2,0],H[2,1],H[2,2])
        str_out += "fractional\n"
        for s in self.structure.symbols:
            for a in self.structure.atomic_basis:
                if a.symbol == s:
                    str_out += "{s} core {r1:10.6f} {r2:10.6f} {r3:10.6f}\n".format(
                            s=s,
                            r1=a.position[0],
                            r2=a.position[1],
                            r3=a.position[2])
        return str_out
    
    def get_gulpinputfile_potentialsection(self):
        str_out = self.potential.gulp_potential_section_to_string(
                parameters=self.parameters)

        return str_out
    
    def restart(self):
        raise NotImplementedError

    def configure_potential(self,configuration=None):
        if configuration is not None:
            self.configuration = copy.deepcopy(configuration)
        #_configuration = self.configuration['potential']
        
        potential_type = self.configuration['potential']['potential_type']
        symbols = self.configuration['potential']['symbols']

        module_name,class_name = potential.PotentialObjectMap(
                potential_type=potential_type)
     
        if potential_type == 'eam':
            eam_setfl_filename = None
            if 'eam_setfl_filename' in _configuraiton:
                self.eam_setfl_filename =_configuration['eam_setfl_filename']
            else:
                # The Embedding Atom model is a special case for the configuration
                # because you do not know the parameterization of the EAM potential
                # unless you know the functional forms of the pair potential,
                # the electron density function, and the embedding function.
                eam_pair_type=_configuration['eam_pair_type']
                eam_density_type=_configuration['eam_density_type']
                eam_embedding_type=_configuration['eam_embedding_type']

                module = importlib.import_module(module_name)
                klass = getattr(module,class_name)
                self.potential = klass(symbols,
                        eam_pair_type,
                        eam_density_type,
                        eam_embedding_type)
        else:
            module = importlib.import_module(module_name)
            klass = getattr(module,class_name)
            self.potential = klass(symbols)
  
        self.parameter_names = list(self.potential.parameter_names)
        self.parameters = self.potential.parameters

    @property
    def parameters(self):
        return self._parameters

    @parameters.setter
    def parameters(self,parameters):
        if self._parameters is None:
            self._parameters = OrderedDict()
        for p in parameters:
            self._parameters[p] = parameters[p]

from pypospack.task.tasks_gulp.gulp_phonons import GulpPhononCalculation

