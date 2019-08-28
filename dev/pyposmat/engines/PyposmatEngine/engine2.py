import os,copy
from collections import OrderedDict
from pypospack.qoi import QoiManager
from pypospack.task import TaskManager
from pypospack.potential import Potential
from pypospack.pyposmat.data.configurationfile import PyposmatConfigurationFile
from pypospack.exceptions import LammpsSimulationError

class PyposmatEngine2(object):
    EAM_EOS_EMBEDDING_FUNCTIONS = ['eam_embed_eos_rose']
    """class for evaluating a simulation

    This class combines the two classes QoiManager, which manages classes used
    to manage the classes used in modelling and calculating material properties,
    and TaskManager, which manages classes used to manage the classes used in
    specific simulations.

    Args:
        filename_in(str): path in
        filename_out(str): path out
        base_directory(str): This is the base directory from which the
            PyposmatEngine will create and run simulations.  By default
            this is set to None, which means it will use the current
            working directory as the base directory.
        fullauto(bool):
    Attributes:
        pyposmat_filename_in(str): path in
        pyposmat_filename_out(str) path out
        base_directory(str): the base directory
        rank_directory(str): This reflect the MPI rank of the processsor
            that the PyposmatEngine is running on.  If there is no MPI
            available, this is automatically set to rank0000.
        configuration(pypospack.pyposmat.PyposmatConfigurationFile)
        qoi_manager(pypospack.qoi.QoiManager)
        task_mamanger(pypospack.task.TaskManager)
    """
    def __init__(self,
            configuration = 'pypospack.config.in',
            data = 'pypospack.results.out',
            base_directory = None,
            fullauto = False):

        self.base_directory = None
        self.configuration = None
        self.qoi_manager = None
        self.task_manager = None

        self._structure_directory = None

        if fullauto:
            self._initialize_configuration(configuration=configuration)
            self._initialize_base_directory(base_directory=base_directory)
            self._initialize_task_manager()
            self._initialize_qoi_manager()

    def _initialize_configuration(self,configuration):
        if instance(configuration,str):
            self.configuration = PyposmatConfigurationFile()
            self.configuration.read(filename=configuration)
        elif isinstance(configuration,PyposmatConfigurationFile):
            self.configuration = configuration
        else:
            msg = ("configuration must either be a string or an instanced "
                   "PyposmatConfigurationFile")
            raise TypeError(msg)

    def _initialize_base_directory(self,base_directory):
        if base_directory is None:
            self.base_directory = os.getcwd()
        elif isinstance(base_directory,str):
            self.base_directory = base_directory
        else:
            msg = "base_directory has to be a string"
            raise TypeError(msg)

        if not os.path.exists(self.base_directory):
            os.mkdir(self.base_directory)

    def _initialize_qoi_manager(self,qois=None):

        if qois is None:
            qois = self.configuration.qois

        self.qoi_manager = QoiManager(qoi_database=qois,fullauto=True)

    def _initialize_task_manager(self):

        self.task_manager = TaskManager(
                base_directory=self.base_directory)
        self.task_manager.configure(
                tasks=self.qoi_manager.tasks,
                structures = self.structures)

    @property
    def structure_directory(self):
        if type(self._structure_directory) is type(str):
            return self._structure_directory
        if type(self.configuration) is not type(None):
            return self.configuration.structures['structure_directory']
        else:
            return None

    @property
    def structures(self):
        """(collections.OrderedDict)"""
        return self.configuration.structures

    @property
    def potential(self):
        """(collections.OrderedDict)"""
        return self.configuration.potential

    def configure(self):
        """

        When writing a new PypospackEngine this method will likely have
        to be modified
        """
        self.create_base_directories()
        self.read_configuration_file()
        self.configure_qoi_manager()
        self.configure_task_manager()

    def create_base_directories(self,base_directory=None):
        assert isinstance(base_directory,str) or base_directory is None

        # <-------- determine the base directory.
        if base_directory is None:
            if self.base_directory is None:
                self.base_directory = os.getcwd()
        elif isinstance(base_directory,str):
            self.base_directory = base_directory
        else:
            msg_err = "the base directory must be a string"
            raise ValueError(msg_err)

        # <-------- create the base directory if the base directory does
        #           not exist
        if not os.path.exists(self.base_directory):
            os.mkdir(self.base_directory)

        # <-------- the rank directory is determined by the MPI rank
        #           this is not implemented yet
        if self.rank_directory is None:
            _rank_directory = "rank0"
            self.rank_directory = os.path.join(
                    self.base_directory,
                    _rank_directory)


    def write_eam_setfl_file(self,parameters,potential,setfl_fn):

        from pypospack.potential import EamPotential
        is_debug = False

        assert isinstance(parameters,OrderedDict)
        assert isinstance(potential,OrderedDict)
        assert isinstance(setfl_fn,str)

        if is_debug:
            print(parameters)
            print(potential)
            print(setfl_fn)

        p = EamPotential(
                symbols=potential['symbols'],
                func_pair=potential['pair_type'],
                func_density=potential['density_type'],
                func_embedding=potential['embedding_type']
        )
        p.write_setfl_file(
                filename=setfl_fn,
                symbols=potential['symbols'],
                Nr=potential['N_r'],
                rmax=potential['r_max'],
                rcut=potential['r_cut'],
                Nrho=potential['N_rho'],
                rhomax=potential['rho_max'],
                parameters=parameters
        )

    def _initialize_potential(potential):
        if isinstance(potential,dict):
            potential_definition = copy.deepcopy(self.configuration.potential)
        elif isinstance(potential,Potential)

    def evaluate_parameter_set(self,parameters):
        """
        Arguments:

            parameters(OrderedDict): hashtable of parameter values with the key value defined
                by the configuration file or the PyposmatConfigurationFile object.  Only the
                free parameter values need to be defined.  Constrained parameter values will
                be automatically generated
        Returns:
            OrderedDict: hashtable of results with the key value being defined
                by the the either configuration file or the PyposmatConfigurationFile object
        Exceptions:
            LammpsSimulationError
        """
        self._initialize_task_manager()
        parameters_ = copy.deepcopy(parameters)
        potential_ = copy.deepcopy(self.configuration.potential)

        # if the filename is EAM, it is faster to create a single EAM file, and
        # then use that file everywhere
        if _potential['potential_type'] == 'eam':
            setfl_fn = os.path.join(
                    os.getcwd(),
                    '{}.eam.alloy'.format("".join(_potential['symbols']))
            )
            if _potential['embedding_type'] in self.EAM_EOS_EMBEDDING_FUNCTIONS:
                self.write_eam_setfl_file(
                        parameters=_parameters,
                        potential=_potential,
                        setfl_fn=setfl_fn
                )
                assert os.path.isfile(setfl_fn)
                _potential['setfl_filename'] = setfl_fn

        try:
            self.task_manager.evaluate_tasks(
                    parameters=_parameters,
                    potential=_potential)
        except LammpsSimulationError as e:
            str_neighlist_overflow = 'Neighbor list overflow'
            raise
        except:
            print("--- FATAL ERROR ---")
            print("self.configuration.potential:")
            for k,v in self.configuration.potential.items():
                print("\t",k,'=',v)
            print("current_parameter_set:")
            for k,v in _parameters.items():
                print("\t",k,'=',v)
            print("--- END ERROR INFO ---")

            print(type(self.configuration.potential))
            raise
        else:
            # send the results from the task calculations to calculate QOIs
            _task_results = self.task_manager.results
            self.qoi_manager.calculate_qois(
                    task_results=_task_results)

            # populate qoi values
            _qoi_results = OrderedDict()
            for k_qoi,v_qoi in self.qoi_manager.qois.items():
                _qoi_val = v_qoi['qoi_val']
                _qoi_results[k_qoi] = _qoi_val

            # populate errors
            _qoi_errors = OrderedDict()
            for k_qoi,v_qoi in self.qoi_manager.qois.items():
                _qoi_error_name = '{}.{}'.format(k_qoi,'err')
                _qoi_error = v_qoi['qoi_err']
                _qoi_errors[_qoi_error_name] = _qoi_error

            _results = OrderedDict()
            _results['parameters'] = copy.deepcopy(_parameters)
            _results['qois'] = copy.deepcopy(_qoi_results)
            _results['errors'] = copy.deepcopy(_qoi_errors)

        return _results
