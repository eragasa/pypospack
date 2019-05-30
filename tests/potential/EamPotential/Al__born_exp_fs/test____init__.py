import pytest
from pypospack.potential import EamPotential

symbols = ['Al']
func_pair_name = "bornmayer"
func_density_name = "eam_dens_exp"
func_embedding_name = "fs"

expected_parameter_names_pair_potential = []
expected_parameter_names_density_function = []
expected_parameter_names_embedding_function = []

expected_parameter_names = [
        'p_AlAl_phi0', 'p_AlAl_gamma', 'p_AlAl_r0', 
        'd_Al_rho0', 'd_Al_beta', 'd_Al_r0', 
        'e_Al_F0', 'e_Al_p', 'e_Al_q', 'e_Al_F1', 'e_Al_rho0']


print(80*'-')
print("func_pair_name={}".format(func_pair_name))
print("func_density_name={}".format(func_density_name))
print("func_embedding_name={}".format(func_density_name))
print(80*'-')

def test____init__():
    obj_pot = EamPotential(
            symbols=symbols,
            func_pair=func_pair_name,
            func_density=func_density_name,
            func_embedding=func_embedding_name)

    assert type(obj_pot) is EamPotential
    assert obj_pot.potential_type == 'eam' 
    assert type(obj_pot.symbols) is list
    assert len(obj_pot.symbols) == len(symbols)
    for i,v in enumerate(symbols):
        obj_pot.symbols[i] = v

    assert obj_pot.is_charge is False
    assert type(obj_pot.parameter_names) is list
    assert len(obj_pot.parameter_names) == len(expected_parameter_names)
    for i,v in enumerate(expected_parameter_names):
        obj_pot.parameter_names = v


if __name__ == "__main__":
    # CONSTRUCTOR TEST

    pot = EamPotential(symbols=symbols,
           func_pair=func_pair_name, 
           func_density=func_density_name,
           func_embedding=func_embedding_name)

    print('pot.potential_type == {}'.format(\
            pot.potential_type))
    print('pot.symbols == {}'.format(\
            pot.symbols))
    print('pot.parameter_names == {}'.format(\
            pot.parameter_names))
    print('pot.is_charge == {}'.format(\
            pot.is_charge))
