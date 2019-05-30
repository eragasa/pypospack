from pypospack.pyposmat import PyposmatDataFile, PyposmatEngine
# from pypospack.qoi import QoiManager
# from pypospack.task import TaskManager

import yaml
from pypospack.io.filesystem import OrderedDictYAMLLoader

class SerialQoiManager(object):
    def __init__(self):
        self.filename_in = None
        self.configuration = None
        self.qois = None
    
    def read(self,filename):
        assert type(filename) is str
        self.filename_in = filename
        _filename = self.filename_in

        _configuration = None
        with open(_filename,'r') as f:
            try:
                _configuration = yaml.load(f,OrderedDictYAMLLoader)
            except yaml.YAMLError as e:
                raise

        self.configuration = OrderdDict()
        self.configuration['qois'] = _configuration['qois']

    def get_task_list(self):
        self.task_list = OrderedDict()
        return self.task_list

    def calculate_qois(self,task_results):

    def get_qois(self):
        pass

    def get_errors(self,ref_type='DFT'):
        pass

class SerialTaskManager(object):
    def __init__(self):
        self.task_list = None
        self.obj_tasks = None
    
    def configure(self,task_list):
        self.task_list = copy.deepcopy(task_list)

        self.obj_tasks = OrderedDict()
        for task_name,task in self.task_list.items:
            self.obj_task['task_name'] = None


class PyposmatEngine(object):
    def __init__(self,
            filename_in = 'pypospack.config.in',
            filename_out = 'pypospack.results.out'):
        self.pyposmat_filename_in = filename_in
        self.pyposmat_filename_out = filename_out

        self.qoi_manager = SerialQoiManager()
        self.task_manager = SerialTaskManager()

    def configure(self):
        self.configure_qoi_manager()
        self.configure_task_manager()
    
    def configure_qoi_manager(self):
        _filename_in = self.pyposmat_filename_in
        self.qoi_manager.read(filename_in=_filename_in)

    def configure_task_manager(self):
        _task_list = self.qoi_manager.get_task_list()
        self.task_manager.configure(
                task_list = _task_list)

    def evaluate_parameter_set(self,parameters):
        
        _parameters = copy.deepcopy(parameters)
        self.task_manager.evaluate_tasks(
                parameters=_parameters)
        _task_results = self.task_manager.results
        self.qoi_manager.calculate_qois(
                task_results=_task_results)
        _qoi_results = self.qoi_manager.results
        _error_results = self.qoi_manager.results

        _results = OrderedDict()
        _results['parameters'] = copy.deepcopy(_parameters)
        _results['qois'] = copy.deepcopy(_qoi_results)
        _results['errors'] = copy.deepcopy(_error_results)

Ni_qoi = OrderedDict():
Ni_qoi['Ni_fcc.E_coh'] = OrderedDict()
Ni_qoi['Ni_fcc.E_coh']['structures'] = {'structure':'Ni_fcc_unit'}
Ni_qoi['Ni_fcc.a1'] = OrderedDict()
Ni_qoi['Ni_fcc.a1']['structures'] = {'structure':'Ni_fcc_unit'}

Ni_potential = OrderedDict()
Ni_potential['potential_type'] = 'eam'
Ni_potential['setfl_filename'] = None
Ni_potential['pair_type'] = 'morse'
Ni_potential['density_type'] = 'eam_dens_exp']
Ni_potential['embedding_type'] = 'eam_embed_universal'
Ni_potential['N_r'] = 10000
Ni_potential['r_max'] = 10.0
Ni_potential['r_cut'] = 10.0
Ni_potential['N_rho'] = 10000
Ni_potential['rho_max'] = 10.0
Ni_potential['symbols'] = ['Ni']

Ni_eam_parameters = OrderedDict()
Ni_eam_parameters['p_NiNi_D0'] = 0.001114
Ni_eam_parameters['p_NiNi_a'] = 3.429506
Ni_eam_parameters['p_NiNi_r0'] = 2.6813
Ni_eam_parameters['d_Ni_rho0'] = 10.0
Ni_eam_parameters['d_Ni_beta'] = 5.0
Ni_eam_parameters['d_Ni_r0'] = 2.0
Ni_eam_parameters['e_Ni_F0'] = 4.10341782e-3
Ni_eam_parameters['e_Ni_p'] = 8.96274624
Ni_eam_parameters['e_Ni_q'] = 8.95940869
Ni_eam_parameters['e_Ni_F1'] = -3.09

_structure_directory = 'test_PypospackEngine'
Ni_structures = OrderedDict()
Ni_structures['Ni_fcc_unit'] = os.path.join(
        _structure_directory,
        'Ni_fcc_unit.gga.relax.vasp')

Ni_pypospack_config = OrderedDict()
Ni_pypospack_config['qoi_info'] = Ni_qoi
Ni_pypospack_config['structures'] = Ni_structures

engine = PyposmatEngine(
        filename_in = 'pypospack.config.in',
        filename_out = 'pypospack.config.out')
engine.configure()

