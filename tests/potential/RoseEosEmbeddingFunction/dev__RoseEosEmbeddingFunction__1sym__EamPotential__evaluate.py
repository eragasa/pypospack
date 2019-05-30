"""
This is a bootstrap script for solving the embedding function implicitly from the parameterization of the embedding and electron density function and
having it being self-consistent with an equation of state

Some of this code was originally developed by CJ OBrien at SANDIA National Laboratory
"""

from pypospack.potential import EamPotential
from collections import OrderedDict
import numpy as np

func_pair_name = "bornmayer"
func_density_name = "eam_dens_exp"
func_embedding_name = "eam_embed_eos_rose"
symbols = ['Ni']

# THE LATTICE INFORMATION IS NORMALLY A QOI, BUT THE QOI'S ARE BURNED HERE
lattice_info = OrderedDict()
for s in symbols:
    lattice_info[s] = OrderedDict()

lattice_info['Ni']['lattice_type'] = 'fcc'
lattice_info['Ni']['cohesive_energy'] = -4.5
lattice_info['Ni']['bulk_modulus'] = 162      # in_GPa
lattice_info['Ni']['lattice_parameter'] = 3.52
a0 = lattice_info['Ni']['lattice_parameter']

# THIS IS COMPUTED INFORMATION AND IS ONLY TRUE FOR AN FCC LATTICE
lattice_type = lattice_info['Ni']['lattice_type']
if lattice_type == 'fcc':

    V = a0**3
    lattice_info['Ni']['equilibrium_volume_per_atom'] = V
    
    re = 1/(2**0.5)*a0
    lattice_info['Ni']['equilibrium_interatomic_distance'] = 1/(2**0.5)*a0 

# PARAMETERS
parameters = OrderedDict()
parameters['p_NiNi_phi0'] = 1.0
parameters['p_NiNi_gamma'] = 2.0
parameters['p_NiNi_r0'] = lattice_info['Ni']['equilibrium_interatomic_distance']
parameters['d_Ni_rho0'] = 1.0
parameters['d_Ni_beta'] = 4.0
parameters['d_Ni_r0'] = lattice_info['Ni']['equilibrium_interatomic_distance']
parameters['e_Ni_ecoh'] = lattice_info['Ni']['cohesive_energy']
parameters['e_Ni_latticetype'] = lattice_info['Ni']['lattice_type']
parameters['e_Ni_B'] = lattice_info['Ni']['bulk_modulus']
parameters['e_Ni_a0'] = lattice_info['Ni']['lattice_parameter']

print(80*'-')
print("func_pair_name={}".format(func_pair_name))
print("func_density_name={}".format(func_density_name))
print("func_embedding_name={}".format(func_embedding_name))
print(80*'-')


# create interatomic distance vector
r_max = 10.
N_r = 1000
r = r_max*np.linspace(1,N_r,N_r)/N_r
print('type(r):{}'.format(type(r)))
print('r:')
print(r)
print('N_r:',r.shape)

# create electron density vector
rho_max = 100000.
N_rho = 1000
rho = rho_max*np.linspace(1,N_rho,N_rho)/N_rho
print('type(rho):{}'.format(type(rho)))
print('rho:')
print(rho)
print('N_rho:',rho.shape)


# CONSTRUCTOR TEST FOR THE EMBEDDING FUNCTION
print(80*'-')
print("CONSTRUCTOR TEST FOR THE EMBEDDING FUNCTION")
print(80*'-')

pot = EamPotential(symbols=symbols,
       func_pair=func_pair_name, 
       func_density=func_density_name,
       func_embedding=func_embedding_name)

from pypospack.potential import EamEmbeddingFunction
from pypospack.potential import EamEmbeddingEquationOfState
assert isinstance(pot.obj_embedding,EamEmbeddingFunction)
assert isinstance(pot.obj_embedding,EamEmbeddingEquationOfState)

pot = EamPotential(symbols=symbols,
       func_pair=func_pair_name, 
       func_density=func_density_name,
       func_embedding=func_embedding_name)

print('type(pot.obj_pair):',type(pot.obj_pair))
print('type(pot.obj_density):',type(pot.obj_density))
print('type(pot.obj_embedding):',type(pot.obj_embedding))
print('pot.obj_pair.parameter_names:',
        pot.obj_pair.parameter_names)
print('pot_obj_density.parameter_names:',
        pot.obj_density.parameter_names)
print('pot.obj_embedding.parameter_names:',
        pot.obj_embedding.parameter_names)

pair_parameters = OrderedDict()
for k,v in parameters.items():
    if k.startswith('p_'):
        pair_parameters[k[2:]] = v
pot.obj_pair.evaluate(
        r=r,
        parameters=pair_parameters
)
print('type(pot.pair):{}'.format(type(pot.pair)))
#for s1 in symbols:
#    for s2 in symbols:
#        k = '{}{}'.format(s1,s2)
#        print('{}:{}:{}'.format(
#            k,
#            type(pot.pair[k]),
#            pot.pair[k].shape))

density_parameters = OrderedDict()
for k,v in parameters.items():
    if k.startswith('d_'):
        density_parameters[k[2:]] = v
pot.obj_density.evaluate(
        r=r,
        parameters=density_parameters
)
#print('type(pot.density):{}'.format(type(pot.density)))
#for s in symbols:
#    print('{}:{}:{}'.format(
#        s,
#        type(pot.density[s]),
#        pot.density[s].shape))


#pot.obj_embedding.evaluate(
#        rho=rho,
#        r=r,
#        parameters=parameters,
#        o_pair=pot.obj_pair,
#        o_density=pot.obj_density)

# RETEST USING AGGREGATED METHOD
print(80*'-')
print("CONSTRUCTOR TEST FOR THE EMBEDDING FUNCTION IN EAMPOTENTIAL CLASS")
print(80*'-')
pot = EamPotential(symbols=symbols,
       func_pair=func_pair_name, 
       func_density=func_density_name,
       func_embedding=func_embedding_name)
pot.evaluate(
        r = r,
        rho = rho,
        rcut = r_max,
        parameters = parameters
    )

print('type(pot.pair):{}'.format(type(pot.pair)))
for s1 in symbols:
    for s2 in symbols:
        k = '{}{}'.format(s1,s2)
        print('{}:{}:{}'.format(
            k,
            type(pot.pair[k]),
            pot.pair[k].shape))

print('type(pot.density):{}'.format(type(pot.density)))
for s in symbols:
    print('{}:{}:{}'.format(
        s,
        type(pot.density[s]),
        pot.density[s].shape))

print('type(pot.embedding):{}'.format(type(pot.embedding)))
for s in symbols:
    print('{}:{}:{}'.format(
        s,
        type(pot.embedding[s]),
      pot.embedding[s].shape))

print('printing these functions')

import matplotlib.pyplot as plt

nearest_neighbor = OrderedDict()
nearest_neighbor[1] = OrderedDict()
nearest_neighbor[2] = OrderedDict()
nearest_neighbor[3] = OrderedDict()
nearest_neighbor[1]['N'] = 12
nearest_neighbor[1]['d/a'] = 1.*np.sqrt(1/2)
nearest_neighbor[2]['N'] = 6
nearest_neighbor[2]['d/a'] = 1.
nearest_neighbor[3]['N'] = 24
nearest_neighbor[3]['d/a'] = 1.*np.sqrt(3/2)

e_coh = 0
e_dens = 0
for i in [1,2,3]:
    a0 = lattice_info['Ni']['lattice_parameter']
    NN_N = nearest_neighbor[i]['N']
    NN_r = nearest_neighbor[i]['d/a'] * a0

    e_coh += np.interp(NN_r,r,pot.pair['NiNi'])
    e_dens += np.interp(NN_r,r,pot.density['Ni'])

e_embed = np.interp(e_dens,rho,pot.embedding['Ni'])

e_coh = e_coh/2 + e_embed
print('E_cohesive={}'.format(e_coh))

fig, ax = plt.subplots(3,1)
for k,v in pot.pair.items():
    ax[0].plot(r,v)
for k,v in pot.density.items():
    ax[1].plot(r,v)
for k,v in pot.embedding.items():
    ax[2].plot(rho,v)

plt.tight_layout()
plt.show()
