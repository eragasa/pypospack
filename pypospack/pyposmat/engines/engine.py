import os,copy
from collections import OrderedDict
from pypospack.qoi import QoiManager
from pypospack.task import TaskManager
from pypospack.pyposmat.data.configurationfile import PyposmatConfigurationFile
from pypospack.task.lammps import LammpsSimulationError

class PyposmatEngine(object):
    """
        Args:
            filename_in(str):
            filename_out(str):
            base_directory(str): This is the base directory from which the
                PyposmatEngine will create and run simulations.  By default
                this is set to None, which means it will use the current
                working directory as the base directory.
            fullauto(bool):
        Attributes:
            pyposmat_filename_in(str)
            pyposmat_filename_out(str)
            base_directory(str)
            rank_directory(str): This reflect the MPI rank of the processsor
                that the PyposmatEngine is running on.  If there is no MPI
                available, this is automatically set to rank0000.
            configuration(pypospack.pyposmat.PyposmatConfigurationFile)
            qoi_manager(pypospack.qoi.QoiManager)
            task_mamanger(pypospack.task.TaskManager)
    """
    def __init__(self,
            filename_in = 'pypospack.config.in',
            filename_out = 'pypospack.results.out',
            base_directory = None,
            fullauto = False):
        assert isinstance(filename_in,str)
        assert isinstance(filename_out,str)
        self.pyposmat_filename_in = filename_in
        self.pyposmat_filename_out = filename_out

        self.base_directory = None
        self.rank_directory = None
        self.configuration = None
        self.qoi_manager = None
        self.task_manager = None

        if base_directory is None:
            self.base_directory = os.getcwd()
        elif isinstance(base_directory,str):
            self.base_directory = base_directory
        else:
            msg_err = "base_directory has to be a string"
            raise ValueError(msg_err)

        if fullauto:
            self.configure()

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
            os.mkdirs(self.base_directory)

        # <-------- the rank directory is determined by the MPI rank
        #           this is not implemented yet
        if self.rank_directory is None:
            _rank_directory = "rank0"
            self.rank_directory = os.path.join(
                    self.base_directory,
                    _rank_directory)


    def read_configuration_file(self,filename=None):
        assert isinstance(filename,str) or filename is None

        _filename_in = None
        if filename is None:
            _filename_in = self.pyposmat_filename_in
        else:
            _filename_in = filename

        self.configuration = PyposmatConfigurationFile(filename=_filename_in)

    def configure_qoi_manager(self,qois=None):
        if qois is not None:
            _qois = qois
        else:
            _qois= self.configuration.qois
        
        self.qoi_manager = QoiManager(qoi_database=_qois,fullauto=True)
    
    def configure_task_manager(self):
        # <-------- local variables
        _base_directory = self.base_directory
        _tasks = self.qoi_manager.tasks
        _structures = self.structures

        # <-------- configure task manager
        self.task_manager = TaskManager(
                base_directory=_base_directory)
        self.task_manager.configure(
                tasks = _tasks,
                structures = _structures)

    def evaluate_parameter_set(self,parameters):
        self.configure_task_manager()
        _parameters = copy.deepcopy(parameters)
        _potential = copy.deepcopy(self.configuration.potential)
        
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
