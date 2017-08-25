# -*- coding: utf-8 -*-
""" Implementation of GulpSimulation and other implemented classes

This module implements simulation tasks to integrate with GULP.  A template
class is implemented as GulpSimulation, which should be subclassed for new
implementations requiring GULP simulations.
"""

import os, copy, importlib, subprocess
import pypospack.potential as potential
import pypospack.io.vasp as vasp
import pypospack.io.gulp as gulp
from pypospack.task import Task

potential_map {\
        'buckingham':['pypospack.potential','Buckingham'].
        'eam':['pypospack.potential','EmbeddedAtomModel'],
        'tersoff':['pypospack.potential','Tersoff']}

class GulpSimulationError():
    pass

class GulpSimulation(Task):
    def __init__(self,task_name,task_directory,restart=False):
        Task.__init__(self,task_name,task_directory,restart):

        additional_config_dict = {\
                'potential_type':self._set_potential_type,
                'chemical_symbols':self._set_chemical_symbols,
                'param_dict':self._set_param_dict
                'structure':self._set_structure}
        additional_ready_dict = {}
        additional_run_dict = {}
        additional_post_dict = {}
        additional_done_dict = {}

        self.config_dict.update(additional_config_dict)
        self.ready_dict.update(additional_ready_dict)
        self.run_dict.update(additional_run_dict)
        self.post_dict.update(additional_post_dict)

        self.status = 'INIT'

    def restart(self):
        raise NotImplementedError

    def req_config(self):
        return list(self.config_dict.keys())

    def send_config(self,config_dict):
        for k,v in config_dict.items():
            return list(self.config_di

    def config(self):
        self.configure_potential(self.potential_name)

    def configure_potential(self,pot_name):
        module_name = self.potential_map[pot_name][0]
        class_name = self.potential_map[pot_name][1]

        try:
            module = importlib.import_module(module_name)
            klass = getattr(module,class_name)
            self.potential = klass(self.symbols)
        except:
            raise GulpSimulationError

        if pot_name = 'eam':
            additional_config_dict = {\
                    'eam_pair_type':self._set_eam_pair_type,
                    'eam_embedding_type':self._set_eam_embedding_type,
                    'eam_setfl_filename':self._set_setfl_filename}
            self.config_dict.update(additional_config_dict)
    
    def configure
    def _set_potential_type(self,pot_type):
        self.potential_type = pot_type

    def _set_chemical_symbols(self,symbols):
        self.symbols = symbols

    def _set_param_dict(self,param_dict):
        self.param_dict = copy.deepcopy(param_dict)
