from collections import OrderedDict
from pypospack.qoi import Qoi


class StackingFaultEnergyCalculation(Qoi):
    def __init__(self, qoi_name, structures):
        qoi_type = "stacking_fault_energy"
        Qoi.__init__(self.qoi_name,qoi_type,structures)
        self._req_var_names = []
        self._req_var_names.append("{}.{}".format(structures[0],'E_min'))
        self._req_var_names.append("{}.{}".format(structures[0],'n_atoms'))
        self._req_var_names.append("{}.{}".format(structures[1],'a1'))
