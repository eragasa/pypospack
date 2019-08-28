import copy, yaml, importlib
from collections import OrderedDict
import numpy as np

def get_qoi_map():
    qoi_map = {
            'gulp_opti_calc':{
                'qoi':['Ecoh_min_all_g',
                       'a11_min_all_g', 'a12_min_all_g','a13_min_all_g',
                       'a21_min_all_g', 'a22_min_all_g','a23_min_all_g',
                       'a31_min_all_g', 'a32_min_all_g','a33_min_all_g',
                       ],
                'module':'pypospack.qoi',
                'class':'GulpOptiCalculations'},
            'lmps_min_none':{
                'qoi':[
                    'Ecoh_min_none',
                    'a1_min_none','a2_min_none','a3_min_none',
                    'a11_min_none','a12_min_none','a13_min_none',
                    'a21_min_none','a22_min_none','a23_min_none',
                    'a31_min_none','a32_min_none','a33_min_none',
                    'p_11_min_none','p_12_min_none','p_13_min_none',
                    'p_21_min_none','p_22_min_none','p_23_min_none',
                    'p_31_min_none','p_32_min_none','P_33_min_none',
                      ],
                'module':'pypospack.qoi',
                'class':'StaticStructureCalculations'},
            'lmps_min_pos':{
                'qoi':['Ecoh_min_pos',
                       'a1_min_pos', 'a2_min_pos', 'a3_min_pos',
                       'a11_min_pos', 'a12_min_pos','a13_min_pos',
                       'a21_min_pos', 'a22_min_pos','a23_min_pos',
                       'a31_min_pos', 'a32_min_pos','a33_min_pos',
                      ],
                'module':'pypospack.qoi',
                'class':'RelaxedPositionCalculations'},
            'lmps_min_all':{
                'qoi':['Ecoh_min_all',
                       'a1_min_all', 'a2_min_all', 'a3_min_all',
                       'a11_min_all', 'a12_min_all','a13_min_all',
                       'a21_min_all', 'a22_min_all','a23_min_all',
                       'a31_min_all', 'a32_min_all','a33_min_all',
                      ],
                'module':'pypospack.qoi',
                'class':'RelaxedStructureCalculations'},
            'lmps_elastic':{
                'qoi':['c11','c12','c13','c22','c23',
                       'c33','c44','c55','c66',
                       'bulk_modulus',
                       'shear_modulus'],
                'module':'pypospack.qoi',
                'class':'ElasticPropertyCalculations'},
            'lmps_phase_order':{
                'qoi':['phase_order'],
                'module':'pypospack.qoi',
                'class':'PhaseOrderCalculation'},
            'lmps_defect':{
                'qoi':['E_formation'],
                'module':'pypospack.qoi',
                'class':'DefectFormationEnergy'},
            'lmps_surface_energy':{
                'qoi':['E_surface'],
                'module':'pypospack.qoi',
                'class':'SurfaceEnergyCalculation'},
            'lmps_stacking_fault':{
                'qoi':['E_stacking_fault'],
                'module':'pypospack.qoi',
                'class':'StackingFaultEnergyCalculation'},
            'lmps_thermal_expansion':{
                'qoi':['thermal_expansion_coefficient'],
                'module':'pypospack.qoi',
                'class':'ThermalExpansion'},
            'gulp_gamma_phonons':{
                'qoi':['gamma_phonon_{}'.format(i) for i in range(20)],
                'module':'pypospack.qoi',
                'class':'GammaPointPhonons'}
        }
    return copy.deepcopy(qoi_map)

def get_supported_qois():
    qois = []
    for k,v in get_qoi_map().items():
        qois += v['qoi']
    return list(qois)

class Qoi(object):
    """ Abstract Quantity of Interest

    Args:
        qoi_name(str)
        qoi_type(str)
        structure_names(list of str)

    Attributes:
        qoi_name(str): this is a unique identifier
        qoi_type(str): this should be set by the implementing class
        structure_names
        reference_values(dict): these are the reference values of quantities
            of interest which we want to calculate.
        task_definitions(dict):
             task_definition[task_name]
             task_definition[task_name][task_type]
             task_definition[task_name][structure]
        task_dependencies(dict):
             task_dependencies[task_name]
        results(dict): these are the results after aggregating the values
            from the different simulations and making some calculations.
    """
    def __init__(self, qoi_name, qoi_type, structures):
        assert isinstance(qoi_name,str)
        assert isinstance(qoi_type,str)
        assert any([
            isinstance(structures,dict),
            isinstance(structures,list),
            isinstance(structures,str),
            ])
        #assert all([isinstance(k,str) for k,v in structures.items()])
        #assert all([isinstance(v,str) for k,v in structures.items()])

        self.qoi_name = qoi_name
        self.qoi_type = qoi_type
        self.structures = copy.deepcopy(structures)
        self.tasks = None

    def determine_tasks(self):
        raise NotImplementedError

    def calculate_qoi(self):
        raise NotImplementedError

    @property
    def reference_value(self):
        return self._ref_value

    @reference_value.setter
    def reference_value(self, value):
        assert type(value), float
        self._ref_value = value

    @property
    def predicted_value(self):
        return self._predicted_value

    @predicted_value.setter
    def predicted_value(self, qhat):
        self._predicted_value = qhat

    @property
    def task_configurations(self):
        return self.tasks

    def __initialize_task_configuration(self):
        if self.tasks is None:
            self.tasks = OrderedDict()

    def __initialize_task_results(self):
        if self.task_results is None:
            self.task_results = OrderedDict()

    def add_task(self,
            task_type,
            task_name,
            task_structure,
            bulk_structure_name=None,
            task_options=None):

        self.__initialize_task_configuration()

        # add task
        self.tasks[task_name] = OrderedDict()
        self.tasks[task_name]['task_type'] = task_type
        self.tasks[task_name]['task_structure'] = task_structure
        if bulk_structure_name is not None:
            self.tasks[task_name]['bulk_structure'] = bulk_structure_name
        if task_options is not None:
            self.tasks[task_name]['task_options'] = task_options

    def process_task_results(self,task_results):
        assert isinstance(task_results,dict)

        self.__initialize_task_results()

    def add_required_simulation(self,structure,simulation_type):
        simulation_name = '{}.{}'.format(structure,simulation_type)
        self.required_simulations[simulation_name] = {}
        self.required_simulations[simulation_name]['structure'] = structure
        self.required_simulations[simulation_name]['simulation_type'] = simulation_type
        self.required_simulations[simulation_name]['precedent_tasks'] = None
        return simulation_name

    def add_precedent_task(self,
            task_name,precedent_task_name,precedent_variables):
        """add a precedent task

        Args:
            task_name(str):
            precedent_task_name(str):
            precedent_variables(str):

        Raises:
            ValueError

        """

        if task_name not in self.required_simulations.keys():
            s = ( 'Tried to add precedent_task_name to task_name.  task_name '
                  'does not exist in required_simulations.\n'
                  '\ttask_name: {}\n'
                  '\tpredecessor_task_name: {}\n' ).format(
                          task_name,
                          precedent_task_name)
            raise ValueError(s)

        if precedent_task_name not in self.required_simulations.keys():
            s = ( 'Tried to add predecessor_task_name to task_name.  '
                  'predecessor_task_name task_name does not exist in '
                  'required_simulations.\n'
                  '\ttask_name: {}\n'
                  '\tpredecessor_task_name: {}\n' ).format(
                          task_name,
                          precedent_task_name)
            raise ValueError(s)

        if self.required_simulations[task_name]['precedent_tasks'] is None:
            self.required_simulations[task_name]['precedent_tasks'] = {}

        self.required_simulations[task_name]['precedent_tasks'][precedent_task_name] ={}
        for v in required_variables:
            self.required_simulations[task_name]['precedent_tasks'][precedent_task_name] = {
                    'variable_name':v,
                    'variable_value':None}

    def determine_required_simulations(self):
        raise NotImplementedError

    def calculate_qoi(self):
        s_type = str(type(self))
        raise NotImplementedError(s_type)

    def get_required_simulations(self):
        if self.required_simulations is None:
            self.determine_required_simulations()
        return copy.deepcopy(self.required_simulations)

from pypospack.qois.crystalstructuregeometry import RelaxedStructureCalculations
from pypospack.qois.crystalstructuregeometry import RelaxedPositionCalculations
from pypospack.qois.lammps_min_none import StaticStructureCalculations
from pypospack.qois.lammps_elastic_properties import ElasticPropertyCalculations
from pypospack.qois.lammps_phase_order import PhaseOrderCalculation
from pypospack.qois.lammps_surface_energy import SurfaceEnergyCalculation
from pypospack.qois.lammps_stacking_fault \
        import StackingFaultEnergyCalculation
from pypospack.qois.lammps_thermal_expansion import ThermalExpansion
from pypospack.qois.lammps_defect_formation_energy \
        import DefectFormationEnergy
from pypospack.qois.gulp_gamma_phonons import GammaPointPhonons
# -----------------------------------------------------------------------------

from pypospack.qois.qoidb import QoiDatabase
from pypospack.qois.qoi_manager import QoiManager

#------------------------------------------------------------------------------
