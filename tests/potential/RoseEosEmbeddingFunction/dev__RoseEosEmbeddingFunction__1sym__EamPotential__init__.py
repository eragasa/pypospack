from pypospack.potential import EamPotential
from collections import OrderedDict
    
# CONSTRUCTOR TEST FOR THE EMBEDDING FUNCTION
symbols = ['Ni']
lattice_info = OrderedDict()
for s in symbols:
    lattice_info[s] = OrderedDict()

func_pair_name = "morse"
func_density_name = "eam_dens_exp"
func_embedding_name = "eam_embed_eos_rose"

lattice_info['Ni']['lattice_type'] = 'fcc'
lattice_info['Ni']['cohesive_energy'] = -4.5
lattice_info['Ni']['bulk_modulus'] = 162      # in_GPa
lattice_info['Ni']['lattice_parameter'] = 3.52
a0 = lattice_info['Ni']['lattice_parameter']

V = a0**3
lattice_info['Ni']['equilibrium_volume_per_atom'] = V

re = 1/(2**0.5)*a0
lattice_info['Ni']['equilibrium_interatomic_distance'] = 1/(2**0.5)*a0 

print(80*'-')
print("func_pair_name={}".format(func_pair_name))
print("func_density_name={}".format(func_density_name))
print("func_embedding_name={}".format(func_embedding_name))
print(80*'-')

# CONSTRUCTOR TEST

pot = EamPotential(symbols=symbols,
       func_pair=func_pair_name, 
       func_density=func_density_name,
       func_embedding=func_embedding_name)

print('pot.potential_type == {}'.format(\
    pot.potential_type))
print('pot.obj_density == {}'.format(
    type(pot.obj_density)))
print('pot.obj_pair_fn == {}'.format(
    type(pot.obj_pair)))
print('pot.obj_embedding_fn == {}'.format(
    type(pot.obj_embedding)))
print('pot.symbols == {}'.format(\
        pot.symbols))
print('pot.parameter_names == {}'.format(\
        pot.parameter_names))
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

