
class AbstractQoi(object):
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

    def add_task(self,task_type,task_name,task_structure,bulk_structure_name=None):
        if self.tasks is None:
            self.tasks = OrderedDict()

        self.tasks[task_name] = OrderedDict()
        self.tasks[task_name]['task_type'] = task_type
        self.tasks[task_name]['task_structure'] = task_structure
        if bulk_structure_name is not None:
            self.tasks[task_name]['bulk_structure'] = bulk_structure_name
    
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
        print(self.qois)
    def configure__get_QoiMap(self):
        self.obj_QoiMap = get_qoi_map()

    def configure__obj_Qoi(self):
        self.obj_qoi = OrderedDict()
        self.qois = OrderedDict()
        for qoik,qoiv in self.qoidb.qois.items():
            self.qois[qoik] = OrderedDict()
            self.qois[qoik]['qoi_name'] = None
            self.qois[qoik]['qoi_ref'] = qoiv['target']

            for qoimapk,qoimapv in self.obj_QoiMap.items():
                if qoiv['qoi_type'] in qoimapv['qoi']:
                    _structures = None
                    _qoi_simulation_type = qoimapk
                    #! determine qoi name
                    if _qoi_simulation_type == 'lmps_phase_order':
                        _s1 = qoiv['structures']['low']
                        _s2 = qoiv['structures']['high']
                        _qoiname = '{}__{}.{}'.format(
                            _s1,_s2,_qoi_simulation_type)
                    elif _qoi_simulation_type in [
                            'lmps_min_all',
                            'lmps_min_none',
                            'lmps_elastic',
                            'lmps_thermal_expansion']:
                        try:
                            _s = qoiv['structures']['ideal']
                        except KeyError as e:
                            print('Qoi:{}'.format(_qoi_simulation_type))
                            print('qoimapk:{}'.format(qoimapk))
                            print('qoik:{}'.format(qoik))
                            print(qoiv['structures'])
                            raise
                        _qoiname = '{}.{}'.format(_s,_qoi_simulation_type)
                    elif _qoi_simulation_type in [
                            'lmps_defect',
                            'lmps_stacking_fault']:
                        _s = qoiv['structures']['defect']
                        _qoiname = '{}.{}'.format(_s,_qoi_simulation_type)
                    elif _qoi_simulation_type in [
                            'lmps_surface_energy']:
                        _s = qoiv['structures']['slab']
                        _qoiname = '{}.{}'.format(_s,_qoi_simulation_type)
                    elif _qoi_simulation_type in [
                            'gulp_gamma_phonons']:
                        _s = qoiv['structures']['ideal']
                        _qoiname = '{}.{}'.format(_s,_qoi_simulation_type)
                    else:
                        msg_err = 'Unknown qoi_simulation_type: {}'
                        msg_err = msg_err.format(_qoi_simulation_type)
                        raise ValueError(msg_err)

                    _qoitype = qoiv['qoi_type']
                    _module = qoimapv['module']
                    _class = qoimapv['class']
                    _structures = qoiv['structures']
                    self._add_obj_Qoi(
                         qoi_name = _qoiname,
                         module_name = _module,
                         class_name = _class,
                         structures = _structures)
                    self.qois[qoik]['qoi_name'] = '{}.{}'.format(_qoiname,_qoitype)
    

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
            task_results(OrderedDict):
        """

        # calculate the QOIs from the QOI objects
        for n_qoi, o_qoi in self.obj_Qoi.items():
            try:
                o_qoi.calculate_qois(task_results=task_results)
            except TypeError as e:
                msg_err = "Cannot calculate qoi for \'{}\'.".format(n_qoi)
                raise ValueError(msg_err+str(e))
        
        
        for k_qoi,_ in self.qois.items():
            try:
                _qoi_id = self.qois[k_qoi]['qoi_name']
                _obj_qoi_id = _qoi_id[:_qoi_id.rindex('.')]
                _qoi_val = self.obj_Qoi[_obj_qoi_id].qois[_qoi_id]
            except:
                print('all_k_qois:')
                for k,v in self.qois.items():
                    print(k,v)
                print('k_qoi:{}'.format(k_qoi))
                print('v_qoi:{}'.format(self.qois[k_qoi]))
                print('obj_qoi_id:{}'.format(_obj_qoi_id))
                print('qoi_id:{}'.format(_qoi_id))
                print(self.obj_Qoi)
                raise
            _qoi_ref = self.qois[k_qoi]['qoi_ref']
            self.qois[k_qoi]['qoi_val'] = _qoi_val
            self.qois[k_qoi]['qoi_err'] = _qoi_val-_qoi_ref
