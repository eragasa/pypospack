import pypospack.potential as potential

symbols = ['Ni']
pot = potential.EamPotential(symbols=symbols)

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

print(80*'-')
symbols = ['Ni','Al']
pot = potential.EamPotential(symbols=symbols)

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

