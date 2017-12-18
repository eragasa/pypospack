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
            'min_none':{
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
                'class':'SinglePointCalculation'},
            'lmps_min_all':{
                'qoi':['Ecoh_min_all',
                       'a1_min_all', 'a2_min_all', 'a3_min_all',
                       'a11_min_all', 'a12_min_all','a13_min_all',
                       'a21_min_all', 'a22_min_all','a23_min_all',
                       'a31_min_all', 'a32_min_all','a33_min_all',
                      ],
                'module':'pypospack.qoi',
                'class':'RelaxedStructureCalculations'},
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

    def add_task(self,task_type,task_name,task_structure):
        if self.tasks is None:
            self.tasks = OrderedDict()

        self.tasks[task_name] = OrderedDict()
        self.tasks[task_name]['task_type'] = task_type
        self.tasks[task_name]['task_structure'] = task_structure

    def process_task_results(self,task_results):
        assert isinstance(task_results,dict)
        
        if self.task_results is None:
            self.task_results = OrderedDict()

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
from pypospack.qois.crystalstructuregeometry import StaticStructureCalculations

# -----------------------------------------------------------------------------
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
    def __init__(self,qoi_database= None,fullauto=True):
        assert any([
            isinstance(qoi_database,QoiDatabase),
            type(qoi_database) is type(None),
            type(qoi_database) is OrderedDict,
            type(qoi_database) is str,
            ])

        #<--------- initialize attributes
        self.qoidb = None
        self.obj_Qoi = None
        self.obj_QoiMap = None
        self.tasks = None

        if isinstance(qoi_database,QoiDatabase):
            self.__init_QoiDatabase_from_QoiDatabase(
                    obj_QoiDatabase=qoi_database)
        elif isinstance(qoi_database,str):
            self.__init_QoiDatabase_from_file(
                    qoi_database_filename=qoi_database)
        elif isinstance(qoi_database,OrderedDict):
            self.__init_QoiDatabase_from_OrderedDict(
                    qoi_database_OrderedDict=qoi_database)
        elif qoi_database is None:
            self.__init_QoiDatabase_from_None()
        else:
            msg_err = (
                "qoi_database argument must be None,QoiDatabase, or str"
                )
            raise ValueError(msg_err)
  
        if fullauto:
            self.configure()
            self.determine_tasks()

    @property
    def qoi_names(self):
        return self.qoidb.qoi_names

    def configure(self):
        self.configure__get_QoiMap()
        self.configure__obj_Qoi()

    def configure__get_QoiMap(self):
        self.obj_QoiMap = get_qoi_map()

    def configure__obj_Qoi(self):
        self.obj_qoi = OrderedDict()

        for qoik,qoiv in self.qoidb.qois.items():
            for qoimapk,qoimapv in self.obj_QoiMap.items():
                if qoiv['qoi_type'] in qoimapv['qoi']:
                    _structures = None
                    _qoi_simulation_type = qoimapk

                    if isinstance(
                            qoiv['structures'],
                            list):
                        _structure = qoiv['structures'][0]
                    elif isinstance(
                            qoiv['structures'],
                            dict):
                        try:
                            _structure = qoiv['structures']['defect']
                        except KeyError:
                            _structure = qoiv['structures']['ideal']
                    elif isinstance(
                            qoiv['structures'],
                            str):
                        _structure = qoiv['structures']
                    else:
                        msg_err = (
                            "Cannot process the type for 'structures':{}"
                            ).format(str(type(qoiv['structures'])))

                    _qoiname = '{}.{}'.format(_structure,_qoi_simulation_type)
                    _module = qoimapv['module']
                    _class = qoimapv['class']
                    _structures = qoiv['structures']
                    self._add_obj_Qoi(
                         qoi_name = _qoiname,
                         module_name = _module,
                         class_name = _class,
                         structures = _structures)

    def calculate_qois(self,task_results):
        assert isinstance(task_results,OrderedDict)

    def __init_QoiDatabase_from_None(self):
        self.qoidb = QoiDatabase()

    def __init_QoiDatabase_from_OrderedDict(self,orderdict_qoi_database):
        raise NotImplementedError

    def __init_QoiDatabase_from_file(self,qoi_database_filename):
        assert isinstance(qoi_database_filename,str)
        self.qoidb = QoiDatabase()
        self.qoidb.read(filename=qoi_database_filename) 

    def __init_QoiDatabase_from_QoiDatabase(self,obj_QoiDatabase):
        assert isinstance(obj_QoiDatabase,QoiDatabase)
        self.qoidb = copy.deepcopy(obj_QoiDatabase)

    def __init_QoiDatabase_from_OrderedDict(self,qoi_database_OrderedDict):
        assert isinstance(qoi_database_OrderedDict,OrderedDict)
        self.qoidb = QoiDatabase(qoi_database_OrderedDict)

    def _add_obj_Qoi(self,qoi_name,module_name,class_name,structures):
        """

        Args:
            qoi_name(str):
            module_name(str):
            class_name(str):
            structures(:obj:`list` of :obj:`str):

        """
        if self.obj_Qoi is None: self.obj_Qoi = OrderedDict()

        if qoi_name not in self.obj_Qoi.keys():
            try:
                module = importlib.import_module(module_name)
                cls = getattr(module,class_name)

                self.obj_Qoi[qoi_name] = cls(
                        qoi_name = qoi_name,
                        structures = structures)
            except:
                print('qoi_name(',type(qoi_name),'):',qoi_name)
                print('module_name(',type(module_name),'):',module_name)
                print('class_name(',type(class_name),'):',class_name)
                print('structures(',type(structures),'):',structures)
                raise

    def determine_tasks(self):
        self.tasks = OrderedDict()

        for k_qoi, v_qoi in self.obj_Qoi.items():
            v_qoi.determine_tasks()

        for k_qoi,v_qoi in self.obj_Qoi.items():
            for k_task, v_task in v_qoi.tasks.items():
                self.add_task(
                        task_name = k_task,
                        task_dict = v_task)

        return copy.deepcopy(self.tasks)

    def add_task(self,task_name,task_dict):
        assert isinstance(task_name,str)
        assert isinstance(task_dict,dict)

        if task_name not in self.tasks:
            self.tasks[task_name] = copy.deepcopy(task_dict)

    def calculate_qois(self,task_results):
        """

        Args:
            results(dict)
        """
        # calculate the material properties from the Qoi objects
        for n_qoi, o_qoi in self.obj_Qoi.items():
            for sim_name, sim_info in o_qoi.get_required_variables().items():
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

#------------------------------------------------------------------------------
from pypospack.io.filesystem import OrderedDictYAMLLoader 
class QoiDatabase(object):
    """ Qoi Database 
    
        Attributes:
            filename(str): file to read/write to yaml file
    """
        
    def __init__(self,qoidb_OrderedDict=None):
        assert any([
            isinstance(qoidb_OrderedDict,OrderedDict),
            type(qoidb_OrderedDict) in [type(None)],
            ])

        self.filename = 'pypospack.qoi.yaml'
        self.qois = None
        self.qoi_names = None

        if qoidb_OrderedDict is not None:
            self.__init_from_OrderedDict(qoidb=qoidb_OrderedDict)

    def __init_from_OrderedDict(self,qoidb):
        for k_qoi,v_qoi in qoidb.items():
            _qoi_name = k_qoi
            _qoi_type = v_qoi['qoi_type']
            _structures = v_qoi['structures']
            _target = v_qoi['target']

            self.add_qoi(
                    qoi_name=_qoi_name,
                    qoi_type=_qoi_type,
                    structures=_structures,
                    target=_target)
        
    def add_qoi(self,
            qoi_name,
            qoi_type,
            structures,
            target):
        """ add a qoi

        Args:
            name(str): name of the qoi.  Usually <structure>.<qoi>.
            qoi(str): name of the qoi.
            structures(list): list of structures
        """
        assert isinstance(qoi_name,str)
        assert isinstance(qoi_type,str)
        assert any([
            isinstance(structures,str),
            isinstance(structures,list),
            isinstance(structures,dict),
            ])
        assert isinstance(target,float)
        
        _structures = None
        if isinstance(structures,list):
            _structures = OrderedDict()
            _structures['ideal'] = structures[0]
        elif isinstance(structures,OrderedDict):
            _structures = copy.deepcopy(structures)
        else:
            raise ValueError

        #<--------- initialize internal atributes if not already set
        if self.qois is None: self.qois = OrderedDict()
        if self.qoi_names is None: self.qoi_names = []

        #<--------- create the dictionary entry for this qoi
        self.qois[qoi_name] = OrderedDict()
        self.qois[qoi_name]['qoi_type'] = qoi_type
        self.qois[qoi_name]['structures'] = copy.deepcopy(_structures)
        self.qois[qoi_name]['target'] = target
        
        #<--------- let's add the value for qoi names
        self.qoi_names.append(qoi_name)

    
    def read(self,filename=None): 
        """ read qoi configuration from yaml file
        Args:
            fname(str): file to read yaml file from.  If no argument is passed
                then use the filename attribute.  If the filename is set, then
                the filename attribute is also set.
        """
        assert isinstance(filename,str)

        # set the attribute if not none
        if filename is not None:
            self.filename = filename
        
        try:
            with open(self.filename,'r') as f:
                self.qois = yaml.load(f, OrderedDictYAMLLoader)
        except:
            raise

        # <------------------ 
        self.qoi_names = [k for k in self.qois]

    def write(self,filename=None):
        """ write qoi configuration to yaml file

        Args:
            fname(str): file to write yaml from from.  If no argument is passed
               then use the filename attribute.  If the filename is set, then 
               the filename attribute is also set.
        """

        # set the attribute if not none
        if filename is not None:
            self.filename = filename

        # marshall attributes into a dictionary
        _qoidb = copy.deepcopy(self.qois)

        # dump to yaml file
        with open(self.filename,'w') as f:
            yaml.dump(_qoidb,f, default_flow_style=False)
#------------------------------------------------------------------------------

class ElasticTensor(Qoi):
    def __init__(self,qoi_name,structures):
        qoi_type = 'elastic_tensor'
        Qoi.__init__(self,qoi_name,qoi_type,structures)
        #self.determine_required_simulations()
        self.structure = self.structures[0]

        required_variables = ['c11','c12','c13','c22','c23','c33',
                              'c44']

    #def determine_required_simulations(self):
    #    if self.required_simulations is not None:
    #        return
    #    self.required_simulations = {}
    #    structure = self.structures[0]
    #    self.add_required_simulation(structure,'elastic')

class DefectFormationEnergy(Qoi):
    def __init__(self,qoi_name, structures):
        if isinstance(structures,list):
            _structures = OrderedDict()
            _structures['defect'] = structures[0]
            _structures['ideal'] = structures[1]
        qoi_type = 'defect_energy'
        Qoi.__init__(self,qoi_name,qoi_type,structures)
        

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
        #self.determine_required_simulations()

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
