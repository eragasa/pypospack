from collections import OrderedDict
from pypospack.qoi import Qoi


class StackingFaultEnergyCalculation(Qoi):
    qois_calculated = ['E_stacking_fault']
    def __init__(self, qoi_name, structures):
        _qoi_name = qoi_name
        _qoi_type = 'lmps_stacking_fault'

        _structures = OrderedDict()
        _structures['ideal'] = structures['ideal']
        _structures['defect'] = structures['defect']
        
        Qoi.__init__(self,
                qoi_name=_qoi_name,
                qoi_type=_qoi_type,
                structures=_structures)

    def determine_tasks(self):
        _ideal_structure_name = self.structures['ideal']
        _ideal_task_type = 'lmps_min_all'
        _ideal_task_name = '{}.{}'.format(
                _ideal_structure_name,
                _ideal_task_type)
        _bulk_structure_name= None
        self.add_task(
                task_type=_ideal_task_type,
                task_name=_ideal_task_name,
                task_structure=_ideal_structure_name,
                bulk_structure_name=_bulk_structure_name)
        
        _defect_structure_name = self.structures['defect']
        _defect_task_type = 'lmps_min_pos'
        _defect_task_name = '{}.{}'.format(
                _defect_structure_name,
                _defect_task_type)
        _bulk_structure_name=self.structures['defect']
        self.add_task(
                task_type=_defect_task_type,
                task_name=_defect_task_name,
                task_structure=_defect_structure_name,
                bulk_structure_name=_bulk_structure_name)

    def calculate_qoi(self):
        _prefix = '{}.{}'.format(
            self.structures['defect'],
            self.qoi_type)
        s_name_defect = self.structures['defect']
        s_name_bulk = self.structures['ideal']
        
        e_slab = task_results[
            "{}.lmps_min_pos.toten".format(s_name_defect)]
        e_bulk = task_results[
            "{}.lmps_min_all.toten".format(s_name_bulk)]
        n_atoms_slab = task_results[
            "{}.lmps_min_pos.natoms".format(s_name_defect)]
        n_atoms_bulk = task_results[
            "{}.lmps_min_all.natoms".format(s_name_bulk)]
        
        a1 = task_results[
            "{}.a11_min_pos".format(s_name_defect)]
        a2 = task_results[
            "{}.a22_min_pos".format(s_name_defect)]
        e_surf = (e_defect - n_atoms_defect/n_atoms_bulk*e_bulk)/(a1*a2)
        
        self.qois = OrderedDict()
        self.qois['{}.E_stacking_fault'.format(_prefix)] = e_surf
