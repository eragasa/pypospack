from collections import OrderedDict
from pypospack.qoi import Qoi

class StaticStructureCalculations(Qoi):
    qois_calculated = [
            'Ecoh_min_none',
            'a1_min_none','a2_min_none','a3_min_none',
            'a11_min_none','a12_min_none','a13_min_none',
            'a21_min_none','a22_min_none','a23_min_none',
            'a31_min_none','a32_min_none','a33_min_none',
            'p_11_min_none','p_12_min_none','p_13_min_none',
            'p_21_min_none','p_22_min_none','p_23_min_none',
            'p_31_min_none','p_32_min_none','P_33_min_none'
    ]

    def __init__(self,qoi_name,structures):
        assert isinstance(qoi_name,str)
        _qoi_name = qoi_name

        _qoi_type = 'lmps_relax_none'

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
        _task_type = 'lmps_min_none'
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
                'lmps_min_none')

        _e_min_none = task_results['{}.{}'.format(_prefix,'toten')]
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

        _a1 = (_a11**2+_a12**2+_a13**2)**0.5
        _a2 = (_a21**2+_a22**2+_a23**2)**0.5
        _a3 = (_a31**2+_a32**2+_a33**2)**0.5

        self.qois = OrderedDict()
        self.qois['{}.Ecoh_min_none'.format(_prefix)] = _e_min_none / _n_atoms
        self.qois['{}.a1_min_none'.format(_prefix)] = _a11
        self.qois['{}.a2_min_none'.format(_prefix)] = _a12
        self.qois['{}.a3_min_none'.format(_prefix)] = _a13
        self.qois['{}.a11_min_none'.format(_prefix)] = _a11
        self.qois['{}.a12_min_none'.format(_prefix)] = _a12
        self.qois['{}.a13_min_none'.format(_prefix)] = _a13
        self.qois['{}.a21_min_none'.format(_prefix)] = _a21
        self.qois['{}.a22_min_none'.format(_prefix)] = _a22
        self.qois['{}.a23_min_none'.format(_prefix)] = _a23
        self.qois['{}.a31_min_none'.format(_prefix)] = _a31
        self.qois['{}.a32_min_none'.format(_prefix)] = _a32
        self.qois['{}.a33_min_none'.format(_prefix)] = _a33
        self.qois['{}.p_11_min_none'.format(_prefix)] = _p11
        self.qois['{}.p_12_min_none'.format(_prefix)] = _p12
        self.qois['{}.p_13_min_none'.format(_prefix)] = _p13
        self.qois['{}.p_21_min_none'.format(_prefix)] = _p21
        self.qois['{}.p_22_min_none'.format(_prefix)] = _p22
        self.qois['{}.p_23_min_none'.format(_prefix)] = _p23
        self.qois['{}.p_31_min_none'.format(_prefix)] = _p31
        self.qois['{}.p_32_min_none'.format(_prefix)] = _p32
        self.qois['{}.p_33_min_none'.format(_prefix)] = _p33
    def get_required_variables(self):
        return list(self.variables.keys())
