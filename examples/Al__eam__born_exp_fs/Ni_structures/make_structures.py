from pypospack.elements import ELEMENTS
from collections import OrderedDict

structures = OrderedDict()

structures['Ni_fcc_unit'] = OrderedDict()
structures['Ni_fcc_unit']['symbols'] = ['Ni']
structures['Ni_fcc_100_unit']['stoichiometry'] = OrderedDict([('Ni',1)])
structures['Ni_fcc_100_unit']['sg'] = 225
structures['Ni_fcc_100_unit']['cell'] = 'unit'

structures['Ni_fcc_orthog'] = OrderedDict()
structures['Ni_fcc_orthog']['symbols'] = ['Ni']
structures['Ni_fcc_orthog']['stoichiometry'] = OrderedDict([('Ni',1)])
structures['Ni_fcc_100_orthog']['sg'] = 225
structures['Ni_fcc_100_orthog']['cell'] = 'orthogonal'

structures['Ni_fcc_primitive'] = OrderedDict()
structures['Ni_fcc_primitive']['symbols'] = ['Ni']
structures['Ni_fcc_primitive']['stoichiometry'] = OrderedDict([('Ni',1)])
structures['Ni_fcc_100_primitive']['sg'] = 225
structures['Ni_fcc_100_primitive']['cell'] = 'primitive'

structures['Ni_fcc_100_unit'] = OrderedDict()
structures['Ni_fcc_100_unit']['symbols'] = ['Ni']
structures['Ni_fcc_100_unit']['stoichiometry'] = OrderedDict([('Ni',1)])
structures['Ni_fcc_100_unit']['sg'] = 225
structures['Ni_fcc_100_unit']['cell'] = [
    [1,0,0],
    [0,1,0],
    [0,0,1]
]

structures['Ni_fcc_110_unit'] = OrderedDict()
structures['Ni_fcc_110_unit']['symbols'] = ['Ni']
structures['Ni_fcc_110_unit']['stoichiometry'] = OrderedDict([('Ni',1)])
structures['Ni_fcc_110_unit']['sg'] = 225
structures['Ni_fcc_110_unit']['lattice_vectors'] = [
    [1,0,0],
    [0,1,0],
    [0,0,1]
]

structures['Ni_fcc_111_unit'] = OrderedDict()
structures['Ni_fcc_111_unit']['symbols'] = ['Ni']
structures['Ni_fcc_111_unit']['stoichiometry'] = OrderedDict([('Ni',1)])
structures['Ni_fcc_111_unit']['sg'] = 225
structures['Ni_fcc_111_unit']['lattice_vectors'] = [
    [1,0,0],
    [0,1,0]
    [0,0,1]
]

structures['Ni_fcc_111_isf'] = OrderedDict()
structures['Ni_fcc_111_isf']['symbols'] = ['Ni']
structures['Ni_fcc_111_isf']['stoichiometry'] = OrderedDict([('Ni',1)])
structures['Ni_fcc_111_isf']['sg'] = 225
structures['Ni_fcc_111_isf']['lattice_vectors'] = [
    [1,0,0],
    [0,1,0]
    [0,0,1]
]
structures['Ni_fcc_111_isf']['defect_type'] = ['stacking_fault']
structures['Ni_fcc_111_isf']['stacking_sequence'] = ['ABC,ABC,AC,ABC,ABC']


structures['Ni_bcc_100_unit'] = OrderedDict()
structures['Ni_bcc_100_unit']['symbols'] = ['Ni']
structures['Ni_bcc_100_unit']['stoichiometry'] = OrderedDict([('Ni',1)])
structures['Ni_bcc_100_unit']['sg'] =
]

class Structure(object):
    def __init__(self, symbols,stoichiometry,sg):
        self.symbols = list(symbols)
        self.stoichiometriy = OrderedDict(stoichiometry)
        self.sg = sg

    def



class StructureDatabase(object):
