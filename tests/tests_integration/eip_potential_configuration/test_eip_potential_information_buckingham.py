import pypospack.potfit as potfit

if __name__ == '__main__':
    pot_info = potfit.PotentialInformation()
    pot_info.potential_type = 'buckingham'
    pot_info.elements = ['Mg','O']
    pot_info.parameter_names = ['chrg_Mg','chrg_O',
                                'MgMg_A','MgMg_rho','MgMg_C',
                                'MgO_A','MgO_rho', 'MgO_C',
                                'OO_A','OO_rho', 'OO_C']
    pot_info.param_info['chrg_Mg'] = {}
    pot_info.param_info['chrg_Mg']['distribution'] = 'uniform'
    pot_info.param_info['chrg_Mg']['low'] = +1.5
    pot_info.param_info['chrg_Mg']['high'] = +2.5
    
    pot_info.param_info['chrg_O'] = {}
    pot_info.param_info['chrg_O']['equals'] = '-chrg_Mg'

    pot_info.param_info['MgO_A'] = {}
    pot_info.param_info['MgO_A']['distribution'] = 'uniform'
    pot_info.param_info['MgO_A']['low'] = 800.0
    pot_info.param_info['MgO_A']['high'] = 1300.0
    
    pot_info.param_info['MgO_rho'] = {}
    pot_info.param_info['MgO_rho']['distribution'] = 'uniform'
    pot_info.param_info['MgO_rho']['low'] = 0.2900
    pot_info.param_info['MgO_rho']['high'] = 0.3300

    pot_info.param_info['MgO_C'] = {}
    pot_info.param_info['MgO_C']['equals'] = 0.00

    pot_info.param_info['OO_A'] = {}
    pot_info.param_info['OO_A']['distribution'] = 'uniform'
    pot_info.param_info['OO_A']['low'] = 500.00
    pot_info.param_info['OO_A']['high'] = 25000.00

    pot_info.param_info['OO_rho'] = {}
    pot_info.param_info['OO_rho']['distribution'] = 'uniform'
    pot_info.param_info['OO_rho']['low'] = 0.1000
    pot_info.param_info['OO_rho']['high'] = 0.4000

    pot_info.param_info['OO_C'] = {}
    pot_info.param_info['OO_C']['distribution'] = 'uniform'
    pot_info.param_info['OO_C']['low'] = 25.0000
    pot_info.param_info['OO_C']['high'] = 77.0000

    pot_info.param_info['MgMg_A'] = {}
    pot_info.param_info['MgMg_A']['equals'] = 0.0

    pot_info.param_info['MgMg_rho'] = {}
    pot_info.param_info['MgMg_rho']['equals'] = 0.5
    
    pot_info.param_info['MgMg_C'] = {}
    pot_info.param_info['MgMg_C']['equals'] = 0.0

    assert pot_info.elements == ['Mg','O']   
    print('elements:',pot_info.elements)
    pot_info.write(fname='pypospack.potential.yaml')

    # check to see if can read what we wrote
    copy_pot_info = potfit.PotentialInformation()
    copy_pot_info.read(fname='pypospack.potential.yaml')

    pot_info = potfit.PotentialInformation()
    pot_info.read(fname='pypospack.buckingham.yaml')

