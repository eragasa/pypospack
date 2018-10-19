
def get_supported_qois():
    supported_qois = ['a0','a1','a2','a3',
                     'alpha','beta','gamma',
                     'c11','c12','c44',
                     'bulk_modulus',
                     'shear_modulus',
                     'defect_energy',
                     'surface_energy',
                     'stacking_fault_energy',
                     'total_energy']
    return supported_qois

def get_supported_potentials():
    supported_potentials = [\
            'buckingham',
            'eam',
            'tersoff']
    return supported_potentials

