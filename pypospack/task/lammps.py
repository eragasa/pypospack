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
import os, copy, shutil,importlib, subprocess
from collections import OrderedDict
import numpy as np
# import pypospack.potential as potential
import pypospack.io.vasp as vasp
import pypospack.io.lammps as lammps
from pypospack.task import Task
from pypospack.io.eamtools import EamSetflFile
from pypospack.potential import Potential,EamPotential,PotentialObjectMap
from pypospack.potential import StillingerWeberPotential

atom_style_list = ['charge','atomic']

lammps_simulation_map = {\
        'min_all':{
            'module':'pypospack.task.lammps',
            'class':'LammpsStructuralMinimization'},
        'min_pos':{
            'module':'pypospack.task.lammps',
            'class':'LammpsPositionMinimization'},
        'min_none':{
            'module':'pypospack.task.lammps',
            'class':'LammpsStaticCalculations'},
        'elastic':{
            'module':'pypospack.task.lammps',
            'class':'LammpsElasticCalculation'},
        'npt':{
            'module':'pypospack.task.lammps',
            'class':'LammpsNptSimulation'}
        }

from pypospack.exceptions import LammpsSimulationError
from pypospack.task.tasks_lammps.abstract_lammps_task import AbstractLammpsSimulation

class LammpsSimulation(AbstractLammpsSimulation): pass

from pypospack.task.tasks_lammps.lmps_min_none import LammpsStaticCalculations
from pypospack.task.tasks_lammps.lmps_min_pos import LammpsPositionMinimization
from pypospack.task.tasks_lammps.lmps_min_all import LammpsStructuralMinimization
from pypospack.task.tasks_lammps.elastic_calculation import LammpsElasticCalculation
from pypospack.task.tasks_lammps.lammps_npt_simulation import LammpsNptSimulation
from pypospack.task.tasks_lammps.lmps_min_sf import LammpsStackingFaultMinimization
