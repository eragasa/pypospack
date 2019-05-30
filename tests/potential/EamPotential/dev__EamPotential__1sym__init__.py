from pypospack.potential import EamPotential

symbols = ['Ni']
func_pair_name = "morse"
func_density_name = "eam_dens_exp"
func_embedding_name = "eam_embed_universal"

print(80*'-')
print("func_pair_name={}".format(func_pair_name))
print("func_density_name={}".format(func_density_name))
print("func_embedding_name={}".format(func_density_name))
print(80*'-')

# CONSTRUCTOR TEST

pot = EamPotential(symbols=symbols,
       func_pair=func_pair_name, 
       func_density=func_density_name,
       func_embedding=func_embedding_name)


print('pot.potential_type == {}'.format(\
        pot.potential_type))
print('pot.symbols == {}'.format(\
        pot.symbols))
print('pot.param_names == {}'.format(\
        pot.param_names))
print('pot.is_charge == {}'.format(\
        pot.is_charge))
print('pot.param == {}'.format(\
        pot.param))
print('type(pot.param) == {}'.format(\
        str(type(pot.param))))

assert pot.potential_type == 'eam'
assert pot.symbols == symbols
assert pot.is_charge is False
assert isinstance(pot.param, dict)
