import copy, yaml
import numpy as np

import copy, importlib

def get_qoi_map():
    qoi_map = {
            'energy':{
                'qoi':['E_coh_sp',
                       'p_xx_sp','p_xy_sp','p_xz_sp',
                       'p_yx_sp','p_yy_sp','p_yz_sp',
                       'p_zx_sp','p_zy_sp','p_zz_sp'],
                'module':'pypospack.qoi',
                'class':'SinglePointCalculation'},
            'crystal':{
                'qoi':['E_coh',
                       'a11','a12','a13',
                       'a21','a22','a23',
                       'a31','a32','a33'],
                'module':'pypospack.qoi',
                'class':'CrystalStructureGeometry'},
            'elastic':{
                'qoi':['c11','c12','c13','c22','c23',
                       'c33','c44','c55','c66',
                       'bulk_modulus',
                       'shear_modulus'],
                'module':'pypospack.qoi',
                'class':'ElasticTensor'},
            'defect_energy':{
                'qoi':['defect_energy'],
                'module':'pypospack.qoi',
                'class':'DefectFormationEnergy'},
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
    qois = []
    for k,v in get_qoi_map().items():
        qois += v['qoi']
    return list(qois)

class Qoi:
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
    def __init__(self, qoi_name, qoi_type, structure_names):
        self.qoi_name = qoi_name
        self.qoi_type = qoi_type
        self.structure_names = list(structure_names)
        
        self.reference_values = None
        self.task_definitions = None
        self.task_dependencies = None
        self.results = None

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
        s_type = str(type(self))
        raise NotImplementedError(s_type)

    def get_required_simulations(self):
        if self.required_simulations is None:
            self.determine_required_simulations()
        return copy.deepcopy(self.required_simulations)

from pypospack.qois.crystalstructuregeometry import CrystalStructureGeometry

class QoiManager(object):
    """ Manager of Quantities of Interest 
    
    This class manages quantities of interest the simulation of multiple 
    quantities of interest will often require the results from the same
    simulation.  The purpose of this class is to identify the values required
    from different simulations, and then identify which simulations need to
    be done first.

    Args:
        qoi_info (str,QoiDatabase,None):  If set to str, then this arguement 
            is treated as a filename used to configure a QoiDatabase object.  If
            QoiDatabase is passed, the reference is set to the qoi_info attribute.
            Default value is None, where the attribute is left as None.
    Attributes:
        qoi_info(pypospack.qoi.QoiDatabase)
        qoi_names (list): list of qoi_names
        qoi_targets (dict): key is the qoi_name, value is the reference value
        qoi_errors (dict): key is the qoi_name, value is the error
        qoi_dict (dict): 
    """ 
    def __init__(self,qoi_info = None):
        self.qoi_map = get_qoi_map()
        self.supported_qois = get_supported_qois()

        self.qoi_info = None
        self.qoi_names = []
        self.obj_qois = {}
         
        if isinstance(qoi_info,QoiDatabase):
            self.qoi_info = qoi_info
            self._config_from_qoi_info()
        elif qoi_info is None:
            pass
        elif isinstance(qoi_info,str):
            self.qoi_info = QoiDatabase(qoi_info)
        else:
            print('qoi_info:',qoi_info,type(qoi_info))
            raise ValueError('qoi_info must be None,QoiDatabase,or string')

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

    def calculate_qois(self,results):
        """

        Args:
            results(dict)
        """
        # calculate the material properties from the Qoi objects
        for qoi, obj_qoi in self.obj_qois.items():
            for sim_name, sim_info in obj_qoi.get_required_variables().items():
                print(qoi,sim_name,sim_info)
        for qoi_name in self.qoi_names:
            print(qoi_name)
            print(self.qoi_info.qois[qoi_name])
        print(80*'-')
        print('required simulations')
        for k,v, in self.required_simulations.items():
            print(k,v)
        for k,v in self.qoi_info.qois.items():
            print(k,v)
        for k,v in self.obj_qois.items():

            print(k,v)

class QoiDatabase(object):
    """ Qoi Database 
   
        Contains methods for managing quantities of interests for a variety 
        of tasks such as marshalling/unmarshalling objects to and from a
        yaml file.  Also includes sanity tests.

        Args:
            filename(str): If this variable is set, this class will attempt to
                unmarhsall the contents of the file into the class.  If this 
                variable is set to None, then the class wil be initialized.
                Default is None.
        Attributes:
            filename(str): file to read/write to yaml file.  By default this
                attribute is set to 'pypospack.qoi.yaml'.
            qois(dict): key is qoi name, value contains a variety of values
    """
        
    def __init__(self, filename = None):
        # initialize attributes
        self.filename = 'pypospack.qoi.yaml'
        self.qois = {}

        if filename is not None:
            self.read(filename)

    @property
    def qoi_names(self):
        return list(self.qois.keys())

    def add_qoi(self,name,qoi,structures,target):
        """ add a qoi

        Args:
            name(str): name of the qoi.  Usually <structure>.<qoi>.
            qoi(str): name of the qoi.
            structures(list): list of structures
        """

        # make a copy of structures
        _structures = None
        if isinstance(structures,str):
            _structures = [structures]
        else:
            _structures = list(structures)

        self.qois[name] = {\
                'qoi':qoi,
                'structures':copy.deepcopy(structures),
                'target':target}

    def read(self,fname=None): 
        """ read qoi configuration from yaml file

        Args:
            fname(str): file to read yaml file from.  If no argument is passed
                then use the filename attribute.  If the filename is set, then
                the filename attribute is also set.
        """

        # set the attribute if not none
        if fname is not None:
            self.filename = fname

        try:
            self.qois = yaml.load(open(self.filename))
        except:
            raise

    def write(self,fname=None):
        """ write qoi configuration to yaml file

        Args:
            fname(str): file to write yaml from from.  If no argument is passed
               then use the filename attribute.  If the filename is set, then 
               the filename attribute is also set.
        """

        # set the attribute if not none
        if fname is not None:
            self.filename = fname

        # marshall attributes into a dictionary
        self.qoi_db = copy.deepcopy(self.qois)

        # dump to yaml file
        with open(self.filename,'w') as f:
            yaml.dump(self.qoi_db,f, default_flow_style=False)

    def check(self):
        """ perform sanity check on potential configuration
        
        does the following checks:
        1.    check to see if the potential type is supported.

        Raises:
            ValueError: if the potential type is unsupported
        # initialize variable, if a check fails the variable will be set to 
        # false
        """

        passed_all_checks = True

        for qoi_name,v in self.qois.items():
            qoi_type = v['qoi']
            if qoi_type not in get_supported_qois():
                print("qoi_name:{}".format(qoi_name))
                print("qoi_type:{}".format(qoi_type))
                print("supported_qois:",get_supported_qois())
                raise ValueError(\
                    "unsupported qoi: {}:{}".format(
                        qoi_name,qoi_type))

    def get_required_structures(self):
        """ get required structures """

        required_structures = []
        for qoi_name, qoi_info in self.qois.items():
            structures = qoi_info['structures']
            for s in structures:
                if s not in required_structures:
                    required_structures.append(s)

        return required_structures


class ElasticTensor(Qoi):
    def __init__(self,qoi_name,structures):
        qoi_type = 'elastic_tensor'
        Qoi.__init__(self,qoi_name,qoi_type,structures)
        self.determine_required_simulations()
        self.structure = self.structures[0]

        required_variables = ['c11','sc12','c13','c22','c23','c33',
                              'c44']
        self.required_variables = {}
        for v in required_variables:
            var_name = "{}.{}.{}".format(
                    self.structure,
                    self.qoi_type,
                    v)
            self.required_variables[var_name] = None

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
        
    def calculate_qoi(self,variables):
        s_name_defect = self.defect_structure
        s_name_bulk   = self.bulk_structure

        #e_defect = self._req_vars["{}.E_min_pos".format(s_name_defect)]
        #e_bulk   = self._req_vars["{}.E_min".format(s_name_bulk)]
        #n_atoms_defect = self._req_vars["{}.n_atoms".format(s_name_defect)]
        #n_atoms_bulk   = self._req_vars["{}.n_atoms".format(s_name_bulk)]
        #e_f = e_defect - n_atoms_defect/n_atoms_bulk*e_bulk
        #self._predicted_value = e_f
        #return self._predicted_value

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
