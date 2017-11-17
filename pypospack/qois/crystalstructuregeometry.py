from collections import OrderedDict
from pypospack.qoi import Qoi

class CrystalStructureGeometry(Qoi):
    def __init__(self,qoi_name,structures):
        self.structure_filenames = OrderedDict()
        qoi_type = 'crystal_structure'
        Qoi.__init__(self,qoi_name,qoi_type,structures)
        self.determine_required_simulations()
        self.variables = {}
        self.variables['xx'] = None
        self.variables['yy'] = None
        self.variables['zz'] = None
        self.variables['xy'] = None
        self.variables['xz'] = None
        self.variables['yz'] = None

    def determine_required_simulations(self):
        if self.required_simulations is not None:
            return

        self.required_simulations = {}
        structure = self.structures[0]
        self.add_required_simulation(structure,'E_min_all')

    def get_required_variables(self):
        return list(self.variables.keys())
