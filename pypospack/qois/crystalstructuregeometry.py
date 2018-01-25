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

    def calculate_qois(self,task_results):
        _prefix = '{}.{}'.format(
                self.structures['ideal'],
                'lmps_min_all')

        _e_min_all = task_results['{}.{}'.format(_prefix,'toten')]
        _n_atoms = task_results['{}.{}'.format(_prefix,'natoms')]
        _p_tot = task_results['{}.{}'.format(_prefix,'totpress')]
        _a11 = task_results['{}.{}'.format(_prefix,'a11')]
        _a12 = task_results['{}.{}'.format(_prefix,'a12')]
        _a13 = task_results['{}.{}'.format(_prefix,'a13')]
        _a21 = task_results['{}.{}'.format(_prefix,'a21')]
        _a22 = task_results['{}.{}'.format(_prefix,'a22')]
        _a23 = task_results['{}.{}'.format(_prefix,'a23')]
        _a31 = task_results['{}.{}'.format(_prefix,'a31')]
        _a32 = task_results['{}.{}'.format(_prefix,'a32')]
        _a33 = task_results['{}.{}'.format(_prefix,'a33')]
        _p11 = task_results['{}.{}'.format(_prefix,'p11')]
        _p12 = task_results['{}.{}'.format(_prefix,'p12')]
        _p13 = task_results['{}.{}'.format(_prefix,'p13')]
        _p21 = task_results['{}.{}'.format(_prefix,'p21')]
        _p22 = task_results['{}.{}'.format(_prefix,'p22')]
        _p23 = task_results['{}.{}'.format(_prefix,'p23')]
        _p31 = task_results['{}.{}'.format(_prefix,'p31')]
        _p32 = task_results['{}.{}'.format(_prefix,'p32')]
        _p33 = task_results['{}.{}'.format(_prefix,'p33')]

        self.qois = OrderedDict()
        self.qois['{}.Ecoh_min_all'.format(_prefix)] = _e_min_all / _n_atoms
        self.qois['{}.a11_min_all'.format(_prefix)] = _a11
        self.qois['{}.a12_min_all'.format(_prefix)] = _a12
        self.qois['{}.a13_min_all'.format(_prefix)] = _a13
        self.qois['{}.a21_min_all'.format(_prefix)] = _a21
        self.qois['{}.a22_min_all'.format(_prefix)] = _a22
        self.qois['{}.a23_min_all'.format(_prefix)] = _a23
        self.qois['{}.a31_min_all'.format(_prefix)] = _a31
        self.qois['{}.a32_min_all'.format(_prefix)] = _a32
        self.qois['{}.a33_min_all'.format(_prefix)] = _a33
    
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
