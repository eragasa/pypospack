from collections import OrderedDict
from pypospack.qoi import Qoi


class StackingFaultEnergyCalculation(Qoi):
    """calculator of stacking fault energies
    
    Args:
       qoi_name(str): the unique identifier of the quantity of interest
       structures(OrderedDict): key is the structure id name, and value is the 
           absolute path of the structure file
    
    """

    qois_calculated = ['E_stacking_fault']
    qoi_type = 'lmps_stacking_fault' 
    def __init__(self, qoi_name, structures):
        assert type(qoi_name) is str
        assert type(structures) is OrderedDict
        assert 'ideal' in structures
        assert 'defect' in structures

        qoi_type = StackingFaultEnergyCalculation.qoi_type

        Qoi.__init__(self,qoi_name=qoi_name,qoi_type=qoi_type,structures=structures)

    def determine_tasks(self):
        """ determine the tasks of the simulation
        
        This method is overridden from the base class.

        Note:
            The first task is structural and atomic position relaxation of 
            ideal bulk.  This task is executed first, and the lattice
            parameters are passed to the second simulation which is
            the calculation of the energy of the stacking fault structure.
        """

        # ideal bulk simulation
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
        
        # stacking fault structure simulation
        _defect_structure_name = self.structures['defect']
        _defect_task_type = 'lmps_min_sf'
        _defect_task_name = '{}.{}'.format(
                _defect_structure_name,
                _defect_task_type)
        _bulk_structure_name=self.structures['ideal']
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
        s_name_bulk = self.structures['ideal']
        
        e_defect = task_results["{}.lmps_min_sf.toten".format(s_name_defect)]
        e_bulk = task_results["{}.lmps_min_all.toten".format(s_name_bulk)]
        n_atoms_defect = task_results["{}.lmps_min_sf.natoms".format(s_name_defect)]
        n_atoms_bulk = task_results["{}.lmps_min_all.natoms".format(s_name_bulk)]
        
        a1 = task_results["{}.lmps_min_sf.a11".format(s_name_defect)]
        a2 = task_results["{}.lmps_min_sf.a22".format(s_name_defect)]
        e_stack = (e_defect - n_atoms_defect/n_atoms_bulk*e_bulk)/(a1*a2)
        
        self.qois = OrderedDict()
        self.qois['{}.E_stacking_fault'.format(_prefix)] = e_stack
