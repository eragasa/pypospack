from collections import OrderedDict
from pypospack.qoi import Qoi

class DefectFormationEnergy(Qoi):

    qois_calculated = ['E_formation']

    def __init__(self,qoi_name,structures):
        assert isinstance(qoi_name,str)
        assert isinstance(structures,dict)

        _qoi_name = qoi_name
        _qoi_type = 'lmps_defect'

        _structures = OrderedDict()
        _structures['ideal'] = structures['ideal']
        _structures['defect'] = structures['defect']

        Qoi.__init__(self,
                qoi_name=_qoi_name,
                qoi_type=_qoi_type,
                structures=_structures)

    def determine_tasks(self):

        # 1. minimize the prototype bulk structure
        _ideal_structure_name = self.structures['ideal']
        _ideal_task_type = 'lmps_min_all'
        _ideal_task_name = '{}.{}'.format(
                _ideal_structure_name,
                _ideal_task_type)
        _bulk_structure_name = None
        self.add_task(
                task_type=_ideal_task_type,
                task_name=_ideal_task_name,
                task_structure=_ideal_structure_name,
                bulk_structure_name=_bulk_structure_name)

        # 2. position minimization of the defect structre
        # '{}.a11_min_all' is used as the length of the a1 vector
        _defect_structure_name = self.structures['defect']
        _defect_task_type = 'lmps_min_pos'
        _defect_task_name = '{}.{}'.format(
                _defect_structure_name,
                _defect_task_type)
        _a0_task_name = '{}.a11_min_all'.format(_ideal_task_name)
        _bulk_structure_name= self.structures['ideal']
        self.add_task(
                task_type=_defect_task_type,
                task_name=_defect_task_name,
                task_structure=_defect_structure_name,
                bulk_structure_name=_bulk_structure_name)

    def calculate_qois(self,task_results):
        _prefix = '{}.{}'.format(
                self.structures['defect'],
                self.qoi_type)
        s_name_defect = self.structures['defect']
        s_name_bulk   = self.structures['ideal']

        print(task_results)
        e_defect = task_results[
                "{}.lmps_min_pos.toten".format(s_name_defect)]
        e_bulk = task_results[
                "{}.lmps_min_all.toten".format(s_name_bulk)]
        n_atoms_defect = task_results[
                "{}.lmps_min_pos.natoms".format(s_name_defect)]
        n_atoms_bulk = task_results[
                "{}.lmps_min_all.natoms".format(s_name_bulk)]
        e_f = e_defect - n_atoms_defect/n_atoms_bulk*e_bulk
        
        self.qois = OrderedDict()
        self.qois['{}.E_formation_defect'.format(_prefix)]= e_f
