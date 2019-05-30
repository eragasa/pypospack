from pypospack.potential import EamPotential
from collections import OrderedDict
class Lattice:

    def __init__(self):
        self.lattice_type = None
        self.a0 = None
        self.is_unit_cell = True

class FccLattice(Lattice):
    
    def __init__(self):
        Lattice.__init__(self)
        self.lattice_type = 'fcc'

    
# CONSTRUCTOR TEST FOR THE EMBEDDING FUNCTION
symbols = ['Ni']
lattice_info = OrderedDict()
for s in symbols:
    lattice_info[s] = OrderedDict()


lattice_info['Ni']['lattice_type'] = 'fcc'
lattice_info['Ni']['cohesive_energy'] = -4.5
lattice_info['Ni']['bulk_modulus'] = 162      # in_GPa
lattice_info['Ni']['lattice_parameter'] = 3.52
a0 = lattice_info['Ni']['lattice_parameter']

V = a0**3
lattice_info['Ni']['equilibrium_volume_per_atom'] = V

re = 1/(2**0.5)*a0
lattice_info['Ni']['equilibrium_interatomic_distance'] = 1/(2**0.5)*a0 

from pypospack.potential import RoseEquationOfStateEmbeddingFunction
o_embed_fn = RoseEquationOfStateEmbeddingFunction(symbols=symbols)
func_pair_name = "morse"
func_density_name = "eam_dens_exp"
func_embedding_name = "eam_embed_eos_rose"

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

#------------------------------------------------------------------------------
parameters = OrderedDict()
parameters['p_NiNi_D0'] = 1.0
parameters['p_NiNi_a'] = 1.0
parameters['p_NiNi_r0'] = lattice_info['Ni']['equilibrium_interatomic_distance']
parameters['d_Ni_rho0'] = 1.0
parameters['d_Ni_beta'] = 1.0
parameters['d_Ni_r0'] = lattice_info['Ni']['equilibrium_interatomic_distance']
parameters['e_Ni_ecoh'] = lattice_info['Ni']['cohesive_energy']
parameters['e_Ni_latticetype'] = lattice_info['Ni']['lattice_type']
parameters['e_Ni_B'] = lattice_info['Ni']['bulk_modulus']
parameters['e_Ni_a0'] = lattice_info['Ni']['lattice_parameter']

pot.evaluate(
        r = np.array(),
        rho = np.array(),
        parameters = parameters
    )

