import copy
import numpy as np

import copy, importlib
import pypospack.potfit as potfit

def get_qoi_map():
    qoi_map = {
            'crystal':{
                'qoi':['a0','a1','a2','a3','alpha','beta','gamma'],
                'module':'pypospack.qoi',
                'class':'CrystalStructureGeometry'},
            'elastic':{
                'qoi':['c11','c12','c13','c22','c23',
                       'c33','c44','c55','c66'],
                'module':'pypospack.qoi',
                'class':'ElasticTensor'},
            'defect_energy':{
                'qoi':['defect_energy'],
                'module':'pypospack.qoi',
                'class':'DefectFormationEnergy'},
            'bulk_modulus':{
                'qoi':['bulk_modulus'],
                'module':'pypospack.qoi',
                'class':'ShearModulus'},
            'surface_energy':{
                'qoi':['surface_energy'],
                'module':'pypospack.qoi',
                'class':'SurfaceEnergy'},
            'stacking_fault_energy':{
                'qoi':['StackingFaultEnergy'],
                'module':'pypospack.qoi',
                'class':'StackingFaultEnergy'}}
    return copy.deepcopy(qoi_map)

def get_supported_qois():
    return list(get_qoi_map().keys())

class QoiManager(object):
    """ Manager of Quantities of Interest 
    
    This class manages quantities of interest the simulation of multiple 
    quantities of interest will often require the results from the same
    simulation.  The purpose of this class is to identify the values required
    from different simulations, and then identify which simulations need to
    be done first.

    Args:

    Attributes:
        qoi_names (list): list of qoi_names
        qoi_targets (dict): key is the qoi_name, value is the reference value
        qoi_errors (dict): key is the qoi_name, value is the error
        qoi_dict (dict): 
    """ 
    def __init__(self,qoi_info = None):
        self.qoi_map = get_qoi_map()
        self.supported_qois = get_supported_qois()

        self.qoi_names = []
        self.obj_qois = {}

        if qoi_info is not None:
            assert isinstance(qoi_info,potfit.QoiDatabase)
            self.qoi_info = qoi_info
            self._config_from_qoi_info()

    def _config_from_qoi_info(self,qoi_info = None):
        """ configure qoi from QoiInformation object

        """

        if qoi_info is not None:
            assert isinstance(potfit.QoiDatabase,qoi_info)
            self.qoi_info = copy.deepcopy(qoi_info)

        self.obj_qois = {}
        self.qoi_names = self.qoi_info.qoi_names

        for qoi, qoi_info in self.qoi_info.qois.items():
            # iterate over items in the qoi configuration dictionary
            for map_qoi, map_info in self.qoi_map.items():
                 # iterate over items in the qoi map dictionary
                 if qoi_info['qoi'] in map_info['qoi']:
                     try:
                         qoi_name = '{}.{}'.format(
                                 qoi_info['structures'][0],
                                 map_qoi)
                     except KeyError as e:
                         print('qoi_info:',qoi_info)
                         raise
                     module_name = map_info['module']
                     class_name = map_info['class']
                     structures = list(qoi_info['structures'])
                     self._add_qoi_to_obj_dict(
                             qoi_name = qoi_name,
                             module_name = module_name,
                             class_name = class_name,
                             structures = structures)

    def _add_qoi_to_obj_dict(self,qoi_name,module_name,class_name,structures):
        """

        Args:
            qoi_name(str):
            module_name(str):
            class_name(str):
            structures(:obj:`list` of :obj:`str):

        """
        if qoi_name not in self.obj_qois.keys():
            try:
                module = importlib.import_module(module_name)
                klass = getattr(module,class_name)

                self.obj_qois[qoi_name] = klass(
                        qoi_name = qoi_name,
                        structures = structures)
            except:
                print('qoi_name(',type(qoi_name),'):',qoi_name)
                print('module_name(',type(module_name),'):',module_name)
                print('class_name(',type(class_name),'):',class_name)
                print('structures(',type(structures),'):',structures)
                raise

    def get_required_simulations(self):
        self.required_simulations = {}
        for k,v in self.obj_qois.items():
            for sim_name,sim_info in v.get_required_simulations().items():
                if sim_name not in self.required_simulations.keys():
                    self.required_simulations[sim_name] = copy.deepcopy(sim_info)
        return copy.deepcopy(self.required_simulations)


class Qoi:
    """ Abstract Quantity of Interest

    Args:
        qoi_name(str)
        qoi_type(str)
        structure_names(list of str)

    Attributes:
        qoi_name(str)
        qoi_type(str)
        ref_value(float)
        required_simulations(dict)
    """
    def __init__(self, qoi_name, qoi_type, structures):
        self.qoi_name = qoi_name
        self.qoi_type = qoi_type

        self.required_simulations = None
        self.predecessor_task = {}
        self.required_variables = {}
        self.ref_value = None
        self.predicted_value = None
        self.error = None

        # set required structure names and set up the dictionary
        self.structures = list(structures)

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
    def set_variable(self, v_name, v_value):
        self._req_vars[v_name] = v_value

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
        raise NotImplementedError

    def get_required_simulations(self):
        if self.required_simulations is None:
            self.determine_required_simulations()
        return copy.deepcopy(self.required_simulations)

class CrystalStructureGeometry(Qoi):
    def __init__(self,qoi_name,structures):
        qoi_type = 'crystal_structure'
        Qoi.__init__(self,qoi_name,qoi_type,structures)
        self.determine_required_simulations()

    def determine_required_simulations(self):
        if self.required_simulations is not None:
            return

        self.required_simulations = {}
        structure = self.structures[0]
        self.add_required_simulation(structure,'E_min_all')


class ElasticTensor(Qoi):
    def __init__(self,qoi_name,structures):
        qoi_type = 'elastic_tensor'
        Qoi.__init__(self,qoi_name,qoi_type,structures)
        self.determine_required_simulations()

    def determine_required_simulations(self):
        if self.required_simulations is not None:
            return
        self.required_simulations = {}
        structure = self.structures[0]
        self.add_required_simulation(structure,'elastic')

class CohesiveEnergy(Qoi):
    def __init__(self,qoi_name,structures):
        qoi_type = 'E_coh'
        Qoi.__init__(self,qoi_name,qoi_type,structures)

    def determine_required_simulations(self):
        if self.required_simulations is not None:
            return
        self.required_simulations = {}
        structure = self.structures[0]
        self.add_required_simulation(structure,'E_min_all')

class BulkModulus(Qoi):
    def __init__(self,qoi_name, structures):
        qoi_type = 'bulk_modulus'
        Qoi.__init__(self,qoi_name,qoi_type,structures)
        self.determine_required_simulations()

    def determine_required_simulations(self):
        if self.required_simulations is not None:
            return
        self.required_simulations = {}
        structure = self.structures[0]
        self.add_required_simulation(structure,'elastic')

class ShearModulus(Qoi):
    def __init__(self,qoi_name, structures):
        qoi_type = 'shear_modulus'
        Qoi.__init__(self,qoi_name,qoi_type,structures)
        self.determine_required_simulations()

    def determine_required_simulations(self):
        if self.required_simulations is not None:
            return
        self.required_simulations = {}
        structure = self.structures[0]
        self.add_required_simulation(structure,'elastic')

class DefectFormationEnergy(Qoi):
    def __init__(self,qoi_name, structures):
        qoi_type = 'defect_energy'
        Qoi.__init__(self,qoi_name,qoi_type,structures)
        self.defect_structure = self.structures[0]
        self.ideal_structure = self.structures[1]
        self.determine_required_simulations()

    def determine_required_simulations(self):
        if self.required_simulations is not None:
            return
        self.required_simulations = {}
        self.add_required_simulation(self.ideal_structure,'E_min_all')
        self.add_required_simulation(self.defect_structure,'E_min_pos')

        # simulation names
        ideal = "{}.{}".format(self.ideal_structure, 'E_min_all')
        defect = "{}.{}".format(self.defect_structure,'E_min_pos')

        precedent_tasks = {}
        precedent_tasks[ideal] = {}
        precedent_tasks[ideal]['variables'] = {
                'a1': None,
                'a2': None,
                'a3': None,
                'E_coh': None,
                'n_atoms': []}

        self.required_simulations[defect]['precedent_tasks'] = \
                copy.deepcopy(precedent_tasks)
        
    def calculate_qoi(self):
        s_name_defect = self._req_structure_names[0]
        s_name_bulk   = self._req_structure_names[1]
        e_defect = self._req_vars["{}.E_min_pos".format(s_name_defect)]
        e_bulk   = self._req_vars["{}.E_min".format(s_name_bulk)]
        n_atoms_defect = self._req_vars["{}.n_atoms".format(s_name_defect)]
        n_atoms_bulk   = self._req_vars["{}.n_atoms".format(s_name_bulk)]
        e_f = e_defect - n_atoms_defect/n_atoms_bulk*e_bulk
        self._predicted_value = e_f
        return self._predicted_value

class SurfaceEnergy(Qoi):
    def __init__(self, qoi_name, structures):
        qoi_type = 'surface_energy'
        Qoi.__init__(self,qoi_name,qoi_type,structures)
        self.surface_structure = self.structures[0]
        self.ideal_structure = self.structures[1]
        self.determine_required_simulations()

    def determine_required_simulations(self):
        if self.required_simulations is not None:
            return

        # adding the simulations
        self.required_simulations = {}
        ideal = self.add_required_simulation(self.ideal_structure,'E_min_all')
        slab = self.add_required_simulation(self.surface_structure,'E_min_pos')

        # define the required simulations
        try:
            self.required_simulations[slab]['precedent_tasks'] = {}
        except KeyError as e:
            s = str(e)
            print(s)
            print('\tslab:{}'.format(slab))
            print('\trequired_simulations:')
            for k,v in self.required_simulations.items():
                print('\t\t',k,v)
        precedent_tasks = {}
        precedent_tasks[ideal] = {}
        precedent_tasks[ideal]['variables'] = {
                'a1': None,
                'a2': None,
                'a3': None,
                'E_coh': None,
                'n_atoms': []}

        self.required_simulations[slab]['precedent_tasks'] = \
                copy.deepcopy(precedent_tasks)

    def calculate_qoi(self):
        s_name_slab = self._req_structure_names[0]
        s_name_bulk = self._req_structure_names[1]
        e_slab = self._req_vars["{}.E_min_pos".format(s_name_slab)]
        e_bulk = self._req_vars["{}.E_min".format(s_name_bulk)]
        n_atoms_slab = self._req_vars["{}.n_atoms".format(s_name_slab)]
        n_atoms_bulk = self._req_vars["{}.n_atoms".format(s_name_bulk)]
        a1 = self._req_vars["{}.a1_min_pos".format(s_name_slab)]
        a2 = self._req_vars["{}.a2_min_pos".format(s_name_slab)]
        e_surf = (e_slab - n_atoms_slab/n_atoms_bulk*e_bulk)/(2*a1*a2)
        self._predicted_value = e_surf
        return self._predicted_value

class StackingFaultEnergy(Qoi):
    def __init__(self, qoi_name, structures):
        qoi_type = "stacking_fault_energy"
        Qoi.__init__(self.qoi_name,qoi_type,structures)
        self._req_var_names = []
        self._req_var_names.append("{}.{}".format(structures[0],'E_min'))
        self._req_var_names.append("{}.{}".format(structures[0],'n_atoms'))
        self._req_var_names.append("{}.{}".format(structures[1],'a1'))
