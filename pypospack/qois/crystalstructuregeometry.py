from collections import OrderedDict
from pypospack.qoi import Qoi

class RelaxedStructureCalculations(Qoi):
    qois_calculated = [
            'a0_min_all',
            'a1_min_all',
            'a2_min_all',
            'a3_min_all',
            'Ecoh_min_all'
    ]
    def __init__(self,qoi_name,structures):
        assert isinstance(qoi_name,str)
        _qoi_name = qoi_name

        _qoi_type = 'lmps_relax_all'

        _structures = OrderedDict()
        if isinstance(structures,str):
            _structures['ideal'] = structures
        elif isinstance(structures,dict):
            _structures['ideal'] = structures['ideal']
        elif isinstance(structures,list):
            _structures['ideal'] = structures[0]
        else:
            msg_err = (
                "structures must be either str, dict, or list"
            )
            raise ValueError(msg_err)

        Qoi.__init__(self,
                qoi_name=qoi_name,
                qoi_type=_qoi_type,
                structures=_structures)

    def determine_tasks(self):
        _structure_ideal_name = self.structures['ideal']
        _task_type = 'lmps_min_all'
        _task_name = "{}.{}".format(
                _structure_ideal_name,
                _task_type)
        _task_requires = None
        self.add_task(
                task_type=_task_type,
                task_name=_task_name,
                task_structure=_structure_ideal_name)
    
    def get_required_variables(self):
        return list(self.variables.keys())

class RelaxedPositionCalculations(Qoi):
    pass
class StaticStructureCalculations(Qoi):
    qois = ['E_min_all.Ecoh',
            'E_min_none.a0',
            'E_min_none.a1',
            'E_min_none.a2',
            'E_min_none.a3',
            'E_min_none.p11',
            'E_min_none.p12',
            'E_min_none.p13',
            'E_min_none.p22',
            'E_min_none,p23',
            'E_min_none,p33']
