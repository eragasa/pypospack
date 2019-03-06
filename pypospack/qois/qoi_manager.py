"""
EJR 09OCT2018 - refactored from pypospack.qoi module, as part of the effort to move all qoi classes as a pypospack.qoi.package, 
"""
import copy,importlib
from collections import OrderedDict
from pypospack.qoi import QoiDatabase
from pypospack.qoi import get_qoi_map

class QoiManager(object):
    """ Manager of Quantities of Interest 
    
    This class manages quantities of interest the simulation of multiple 
    quantities of interest will often require the results from the same
    simulation.  The purpose of this class is to identify the values required
    from different simulations, and then identify which simulations need to
    be done first.

    Args:
        qoi_database(OrderedDict):this is a dictionary of key-value pairs, where 
            the key value is unique qoi_id for this set of simulations, the value 
            contains the necessary configuration information necessary for 
            configuration of object instantiated with Qoi base class.  
            Normally, this object is passed from the PyposmatConfigurationFile
            object using the `qois` attribute.
        fullauto(bool,optional): Set to True by default; this class will be
            initialized using the information in the `qoi_database` argument.
            The purpose of setting this to False is for development purposes
            so that helper methods can be individually tests.
    Attributes:
        qois(OrderedDict):this is a dictionary of key-value pairs, where the key
            value is unique qoi_id for this set of simulations, the value 
            contains the necessary configuration information necessary for 
            configuration of object instantiated with Qoi base class.
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

        self.__init_QoiDatabase(qoi_database=qoi_database)
  
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

    def get_qoi_name(self,qoi_simulation_type,structures):

        if qoi_simulation_type in ['lmps_phase_order']:
            s1 = structures['low']
            s2 = structures['high']
            qoiname = '{}__{}.{}'.format(s1,s2,qoi_simulation_type)

        elif qoi_simulation_type in ['lmps_min_all',
                                     'lmps_min_none',
                                     'lmps_elastic',
                                     'lmps_thermal_expansion',
                                     'gulp_gamma_phonons']:
            s = structures['ideal']
            qoiname = '{}.{}'.format(s,qoi_simulation_type)

        elif qoi_simulation_type in ['lmps_defect',
                                     'lmps_stacking_fault']:
            s = structures['defect']
            qoiname = '{}.{}'.format(s,qoi_simulation_type)

        elif qoi_simulation_type in ['lmps_surface_energy']:
            s = structures['slab']
            qoiname = '{}.{}'.format(s,qoi_simulation_type)

        else:
            msg_err = 'Unknown qoi_simulation_type: {}'.format(
                    qoi_simulation_type)
            raise ValueError(msg_err)
        
        return qoiname
    
    def configure__obj_Qoi(self):
        self.obj_qoi = OrderedDict()
        self.qois = OrderedDict()
        for qoik,qoiv in self.qoidb.qois.items():
            self.qois[qoik] = OrderedDict()
            self.qois[qoik]['qoi_name'] = None
            self.qois[qoik]['qoi_type'] = qoiv['qoi_type']
            self.qois[qoik]['qoi_structures'] = qoiv['structures']
            self.qois[qoik]['qoi_ref'] = qoiv['target']
            try:
                self.qois[qoik]['qoi_options'] = qoiv['qoi_options']
            except KeyError as e:
                pass

            for qoimapk,qoimapv in self.obj_QoiMap.items():
                if qoiv['qoi_type'] in qoimapv['qoi']:
                    _structures = None
                    _qoi_simulation_type = qoimapk
                    _qoiname = self.get_qoi_name(
                        qoi_simulation_type=_qoi_simulation_type,
                        structures=qoiv['structures'])
                    
                    _qoitype = qoiv['qoi_type']
                    _module = qoimapv['module']
                    _class = qoimapv['class']
                    _structures = qoiv['structures']
                    
                    try:
                        _qoi_options = qoiv['qoi_options']
                    except KeyError as e:
                        _qoi_options = None
                    
                    self._add_obj_Qoi(
                         qoi_name = _qoiname,
                         module_name = _module,
                         class_name = _class,
                         structures = _structures,
                         qoi_options = _qoi_options)
                    self.qois[qoik]['qoi_name'] = '{}.{}'.format(_qoiname,_qoitype)

    def __init_QoiDatabase(self,qoi_database=None):

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
    
    def __init_QoiDatabase_from_None(self):
        self.qoidb = QoiDatabase()

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

    def _add_obj_Qoi(self,qoi_name,module_name,class_name,structures,qoi_options=None):
        """

        Args:
            qoi_name(str):
            module_name(str):
            class_name(str):
            structures(:obj:`list` of :obj:`str):

        """
        
        if self.obj_Qoi is None: 
            self.obj_Qoi = OrderedDict()

        if qoi_name not in self.obj_Qoi.keys():
            try:
                module = importlib.import_module(module_name)
                cls = getattr(module,class_name)

                if qoi_options is None:
                    self.obj_Qoi[qoi_name] = cls(
                            qoi_name = qoi_name,
                            structures = structures)
                else:
                    kwargs = qoi_options
                    self.obj_Qoi[qoi_name] = cls(
                            qoi_name = qoi_name,
                            structures = structures,
                            **kwargs)
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
        
        
        for k_qoi,v_qoi in self.qois.items():
            try:
                qoi_id = self.qois[k_qoi]['qoi_name']
                obj_qoi_id = qoi_id[:qoi_id.rindex('.')]
                qoi_val = self.obj_Qoi[obj_qoi_id].qois[qoi_id]
            except:
                print('k_qoi:{}'.format(k_qoi))
                print('qoi_id:{}'.format(qoi_id))
                print('qoi_keys:')
                for k,v in self.qois.items():
                    print("\t{:15}{}".format(k,self.qois[k]['qoi_name']))
                print('all_k_qois:')
                for k,v in self.qois.items():
                    print(k,v)
                print('k_qoi:{}'.format(k_qoi))
                print('v_qoi:{}'.format(self.qois[k_qoi]))
                print('obj_qoi_id:{}'.format(obj_qoi_id))
                print(self.obj_Qoi)
                raise

            qoi_ref = self.qois[k_qoi]['qoi_ref']
            self.qois[k_qoi]['qoi_val'] = qoi_val
            self.qois[k_qoi]['qoi_err'] = qoi_val-qoi_ref
