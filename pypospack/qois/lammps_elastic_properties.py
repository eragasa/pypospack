from collections import OrderedDict
from pypospack.qoi import Qoi

class ElasticPropertyCalculations(Qoi):

    qois_calculated = ['c11','c12','c13','c22','c33','c44','c55','c66']
    
    def __init__(self,qoi_name,structures):
        assert isinstance(qoi_name,str)
        
        _qoi_name = qoi_name
        _qoi_type = 'lmps_elastic'
        
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
        _task_type = 'lmps_elastic'
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
                'lmps_elastic')

        # debugging code
        # if True:
        #    for k,v in task_results.items():
        #        print(k,v)

        _c11 = task_results['{}.{}'.format(_prefix,'c11')]
        _c12 = task_results['{}.{}'.format(_prefix,'c12')]
        _c13 = task_results['{}.{}'.format(_prefix,'c13')]
        _c22 = task_results['{}.{}'.format(_prefix,'c22')]
        _c33 = task_results['{}.{}'.format(_prefix,'c33')]
        _c44 = task_results['{}.{}'.format(_prefix,'c44')]
        _c55 = task_results['{}.{}'.format(_prefix,'c55')]
        _c66 = task_results['{}.{}'.format(_prefix,'c66')]

        # References:
        # http://homepages.engineering.auckland.ac.nz/~pkel015/SolidMechanicsBooks/Part_I/BookSM_Part_I/06_LinearElasticity/06_Linear_Elasticity_03_Anisotropy.pdf

        _bulk = (_c11+2*_c12)/3
        _shear = (_c11-_c12)/2

        self.qois = OrderedDict()
        self.qois['{}.c11'.format(_prefix)] = _c11
        self.qois['{}.c12'.format(_prefix)] = _c12
        self.qois['{}.c13'.format(_prefix)] = _c13
        self.qois['{}.c22'.format(_prefix)] = _c22
        self.qois['{}.c33'.format(_prefix)] = _c33
        self.qois['{}.c44'.format(_prefix)] = _c44
        self.qois['{}.c55'.format(_prefix)] = _c55
        self.qois['{}.c66'.format(_prefix)] = _c66
        self.qois['{}.bulk_modulus'.format(_prefix)] = _bulk
        self.qois['{}.shear_modulus'.format(_prefix)] = _shear
