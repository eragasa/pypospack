from collections import OrderedDict
from pypospack.qoi import Qoi

class SurfaceEnergyCalculation(Qoi):
    def __init__(self, qoi_name, structures):
        qoi_type = 'surface_energy'
        Qoi.__init__(self,qoi_name,qoi_type,structures)
        self.surface_structure = self.structures[0]
        self.ideal_structure = self.structures[1]
        #self.determine_required_simulations()

    def calculate_qoi(self):
        s_name_slab = self._req_structure_names[0]
        s_name_bulk = self._req_structure_names[1]
        e_slab = self._req_vars["{}.E_min_pos".format(s_name_slab)]
        e_bulk = self._req_vars["{}.E_min".format(s_name_bulk)]
        n_atoms_slab = self._req_vars["{}.n_atoms".format(s_name_slab)]
        n_atoms_bulk = self._req_vars["{}.n_atoms".format(s_name_bulk)]
        a1 = self._req_vars["{}.a1_min_pos".format(s_name_slab)]
        a2 = self._req_vars["{}.a2_min_pos".format(s_name_slab)]
        e_surf = (e_slab - n_atoms_slab/n_atoms_bulk*e_bulk)/(2*a1*a2)
        self._predicted_value = e_surf
        return self._predicted_value
