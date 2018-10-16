from collections import OrderedDict
from pypospack.qoi import Qoi

class GammaPointPhonons(Qoi):

    def __init__(self,qoi_name,structures):
        assert isinstance(qoi_name,str)
        assert isinstance(structures,dict)

        _qoi_name = qoi_name
        _qoi_type = 'gulp_gamma_phonons'

        _structures = OrderedDict()
        _structures['ideal'] = structures['ideal']

        Qoi.__init__(self,
                qoi_name=_qoi_name,
                qoi_type=_qoi_type,
                structures=_structures)

    def determine_tasks(self):

        _ideal_structure_name = self.structures['ideal']
        _ideal_task_type = 'gulp_gamma_phonons'
        _ideal_task_name = '{}.{}'.format(
                _ideal_structure_name,
                _ideal_task_type
                )

        _bulk_structure_name = None
        self.add_task(
                task_type=_ideal_task_type,
                task_name=_ideal_task_name,
                task_structure=_ideal_structure_name,
                bulk_structure_name=_bulk_structure_name
                )

    def calculate_qois(self,task_results):
        _prefix = '{}.{}'.format(
                self.structures['ideal'],
                self.qoi_type)

        s_name_ideal = self.structures['ideal']

        self.qois = OrderedDict()
        for k,v in task_results.items():
            print("\t{}:{}".format(k,v))
            self.qois[k] = v

