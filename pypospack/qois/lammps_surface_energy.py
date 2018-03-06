from collections import OrderedDict
from pypospack.qoi import Qoi

class SurfaceEnergyCalculation(Qoi):
    qois_calculated = ['E_surface']
    def __init__(self, qoi_name, structures):
        _qoi_name = qoi_name
        _qoi_type = 'lmps_surface_energy'

        _structures = OrderedDict()
        _structures['ideal'] = structures['ideal']
        _structures['slab'] = structures['slab']
        
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
        
        _slab_structure_name = self.structures['slab']
        _slab_task_type = 'lmps_min_pos'
        _slab_task_name = '{}.{}'.format(
                _slab_structure_name,
                _slab_task_type)
        _bulk_structure_name=self.structures['ideal']
        self.add_task(
                task_type=_slab_task_type,
                task_name=_slab_task_name,
                task_structure=_slab_structure_name,
                bulk_structure_name=_bulk_structure_name)
    
    def calculate_qois(self,task_results):
        _prefix = '{}.{}'.format(
            self.structures['slab'],
            self.qoi_type)
        s_name_slab = self.structures['slab']
        s_name_bulk = self.structures['ideal']
        
        e_slab = task_results[
            "{}.lmps_min_pos.toten".format(s_name_slab)]
        e_bulk = task_results[
            "{}.lmps_min_all.toten".format(s_name_bulk)]
        n_atoms_slab = task_results[
            "{}.lmps_min_pos.natoms".format(s_name_slab)]
        n_atoms_bulk = task_results[
            "{}.lmps_min_all.natoms".format(s_name_bulk)]
        
        a1 = task_results[
            "{}.lmps_min_pos.a11".format(s_name_slab)]
        a2 = task_results[
            "{}.lmps_min_pos.a22".format(s_name_slab)]
        e_surf = (e_slab - n_atoms_slab/n_atoms_bulk*e_bulk)/(2*a1*a2)
        
        self.qois = OrderedDict()
        self.qois['{}.E_surface'.format(_prefix)] = e_surf
