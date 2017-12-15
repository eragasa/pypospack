from collections import OrderedDict

MgO_structures = OrderedDict()
MgO_structures['structure_db_dir'] = 'test_LammpsStructuralMinimization'
MgO_structures['MgO_NaCl_unit'] = OrderedDict()
MgO_structures['MgO_NaCl_unit']['filename'] = 'MgO_NaCl_unit.gga.relax.vasp'
MgO_structures['MgO_NaCl_fr_a'] = OrderedDict()
MgO_structures['MgO_NaCl_fr_a']['filename'] = 'MgO_NaCl_fr_a.gga.relax.vasp'
MgO_structures['MgO_NaCl_fr_c'] = OrderedDict()
MgO_structures['MgO_NaCl_fr_c']['filename'] = 'MgO_NaCl_fr_c.gga.relax.vasp'
MgO_structures['MgO_NaCl_sh'] = OrderedDict()
MgO_structures['MgO_NaCl_sh']['filename'] = 'MgO_NaCl_sh.gga.relax.vasp'
MgO_structures['MgO_NaCl_001s'] = OrderedDict()
MgO_structures['MgO_NaCl_001s']['filename'] = 'MgO_NaCl_001s.gga.relax.vasp'

MgO_LewisCatlow = OrderedDict()
MgO_LewisCatlow['potential'] = OrderedDict()
MgO_LewisCatlow['potential']['potential_type'] = 'buckingham'
MgO_LewisCatlow['potential']['symbols'] = ['Mg','O']

MgO_LewisCatlow['parameters'] = OrderedDict()
# Charge potentials are expected in the order in ['potential']['symbols']
MgO_LewisCatlow['parameters']['chrg_Mg'] = 2.0
MgO_LewisCatlow['parameters']['chrg_O'] = -2.0
# For pair potentials, the order of pairings are determined by the order
# in the ['potential']['symbols'] entry in the dictionary.  In this case,
# MgMg, MgO, and OO
MgO_LewisCatlow['parameters']['MgMg_A'] = 0.0
MgO_LewisCatlow['parameters']['MgMg_rho'] = 0.5
MgO_LewisCatlow['parameters']['MgMg_C'] = 0.0
MgO_LewisCatlow['parameters']['MgO_A'] = 821.6
MgO_LewisCatlow['parameters']['MgO_rho'] = 0.3242
MgO_LewisCatlow['parameters']['MgO_C'] = 0.0
MgO_LewisCatlow['parameters']['OO_A'] = 2274.00
MgO_LewisCatlow['parameters']['OO_rho'] = 0.1490
MgO_LewisCatlow['parameters']['OO_C'] = 27.88


