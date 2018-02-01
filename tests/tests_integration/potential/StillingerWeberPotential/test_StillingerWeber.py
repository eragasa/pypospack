import numpy as np
from collections import OrderedDict

#import pypospack.potential as potential
if __name__ == "__main__":
    symbols = ['Si']
    parameter_names = [
        'SiSiSi_epsilon',
        'SiSiSi_sigma',
        'SiSiSi_a',
        'SiSiSi_lambda',
        'SiSiSi_gamma',
        'SiSiSi_costheta0',
        'SiSiSi_A',
        'SiSiSi_B',
        'SiSiSi_p',
        'SiSiSi_q']

    # Potential parameters from Pizzagalli, J Phys Condens Matter 25 (2013) 055801
    Si_sw_pizzagalli = OrderedDict()
    Si_sw_pizzagalli['SiSiSi_epsilon']=1.04190
    Si_sw_pizzagalli['SiSiSi_sigma']=2.128117
    Si_sw_pizzagalli['SiSiSi_a']=1.80
    Si_sw_pizzagalli['SiSiSi_lambda']=31.0
    Si_sw_pizzagalli['SiSiSi_gamma']=1.10
    Si_sw_pizzagalli['SiSiSi_costheta0']=-1/3
    Si_sw_pizzagalli['SiSiSi_A']=19.0
    Si_sw_pizzagalli['SiSiSi_B']=0.65
    Si_sw_pizzagalli['SiSiSi_p']=3.5
    Si_sw_pizzagalli['SiSiSi_q']=0.5
    
    param_dict = OrderedDict()
    from pypospack.potential import StillingerWeberPotential
    swpot = StillingerWeberPotential(symbols=symbols)
    _sw_parameters = swpot.parameter_names
    swpot.fname_sw_params = ''.join(swpot.symbols) + ".sw.parameters"
    _sw_potential_mod_str = swpot.lammps_potential_section_to_string()
    _sw_potential_param_file_str = swpot.lammps_sw_file_to_string(
            parameters=Si_sw_pizzagalli)
    
    assert type(swpot.potential_type) is str
    assert swpot.potential_type == 'stillingerweber' 
    assert type(swpot.symbols) is list
    assert swpot.symbols == symbols
    assert type(swpot.parameter_names) is list
    assert swpot.parameter_names == parameter_names
    assert type(swpot.is_charge) is bool
    assert swpot.is_charge is False
    
    print("potential_type:{}".format(swpot.potential_type))
    print("symbols:{}".format(swpot.symbols))
    print("param_names:{}".format(swpot.param_names))
    print("is_charge:{}".format(swpot.is_charge))

    print(80*'-')
    print('{:^80}'.format('PARAMETER LIST'))
    print(80*'-')
    for p in _sw_parameters:
        print(p)

    print(80*'-')
    print('{:^80}'.format('POTENTIAL.MOD SECTION'))
    print(80*'-')
    print(_sw_potential_mod_str)
    
    print(80*'-')
    print('{:^80}'.format('SW.POTENTIAL SECTION'))
    print(80*'-')
    print(_sw_potential_param_file_str)

