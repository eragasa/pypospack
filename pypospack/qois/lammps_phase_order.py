from collections import OrderedDict
from pypospack.qoi import Qoi

class PhaseOrderCalculation(Qoi):
    qoi_type = 'lmps_phase_order'
    qois_calculated = ['phase_order']

    def __init__(self,qoi_name,structures):
        assert isinstance(qoi_name,str)
        _qoi_name = qoi_name

        _structures = OrderedDict()
        _structures['low'] = structures['low']
        _structures['high'] = structures['high']

        Qoi.__init__(self,
                qoi_name=qoi_name,
                qoi_type=self.qoi_type,
                structures=_structures)

    def determine_tasks(self):
        _structure_low_name = self.structures['low']
        _task_type = 'lmps_min_all'
        _task_name = '{}.{}'.format(
                _structure_low_name,
                _task_type)
        _task_requires = None
        self.add_task(
                task_type=_task_type,
                task_name=_task_name,
                task_structure=_structure_low_name)

        _structure_high_name = self.structures['high']
        _task_type = 'lmps_min_all'
        _task_name = '{}.{}'.format(
                _structure_high_name,
                _task_type)
        _task_requires = None
        self.add_task(
                task_type=_task_type,
                task_name=_take_name,
                task_structure=_structure_high_name)

    def calculate_qois(self,task_results):
        _low_prefix = '{}.{}'.format(
                self.structures['low'],
                'lmps_min_all')

        _high_prefix = '{}.{}'.format(
                self.structures['high'],
                'lmps_min_all')

        _high_e_min_all = task_results['{}.{}'.format(_high_prefix,'toten')]
        _high_n_atoms = task_results['{}.{}'.format(_high_prefix,'natoms')]
        _high_ecoh = _high_e_min_all/_high_n_atoms

        _low_e_min_all = task_results['{}.{}'.format(_low_prefix,'toten')]
        _low_n_atoms = task_results['{}.{}'.format(_low_prefix,'natoms')]
        _low_ecoh = _low_e_min_all/_low_n_atoms

        _phase_order = _high_ecoh - _low_ecoh

        self.qois = OrderedDict()
        self.qois['{}.{}.phase_order'.format(
            self.structures['low'],
            self.structures['high'])] = _phase_order



