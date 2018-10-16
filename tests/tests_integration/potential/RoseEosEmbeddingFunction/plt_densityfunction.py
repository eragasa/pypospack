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
N_r = 100
r = r_max*np.linspace(1,N_r)/N_r
print('type(r):{}'.format(type(r)))
print('r:')
print(r)
print('N_r:',r.shape)

# create electron density vector
rho_max = 100000.
N_rho = 100
rho = rho_max*np.linspace(1,N_rho)/N_rho
print('type(rho):{}'.format(type(rho)))
print('rho:')
print(rho)
print('N_rho:',rho.shape)

pot = EamPotential(symbols=symbols,
       func_pair=func_pair_name, 
       func_density=func_density_name,
       func_embedding=func_embedding_name)

from pypospack.potential import EamEmbeddingFunction
from pypospack.potential import EamEmbeddingEquationOfState
assert isinstance(pot.obj_embedding,EamEmbeddingFunction)
assert isinstance(pot.obj_embedding,EamEmbeddingEquationOfState)

def rhofxn(a,args):
    rho_fxn_parameters = OrderedDict()
    for k,v in args[0].items():
        if k.startswith('d_{}'.format(s)):
            rho_fxn_parameters[k[2:]] = v

    lattice_type = args[0]['e_{}_latticetype'.format(s)]
    if lattice_type == 'fcc':
        n_NN = [12,6,24,12,24,8]
        d_NN = [a/np.sqrt(2.),a,a*np.sqrt(1.5),a*np.sqrt(2.0),a*np.sqrt(2.5),a*np.sqrt(3.0)]
    
    _rho = 0.
    for i in range(len(n_NN)):
        if len(rho_fxn_parameters) == 0:
            _rho += n_NN[i] * pot.obj_density.evaluate(
                    d_NN[i],
                    rho_fxn_parameters
                )[s]
        else:
            _rho += n_NN[i] * pot.obj_density.evaluate(
                    d_NN[i],
                    rho_fxn_parameters
                )[s]
    return _rho

N_lattice_parameter = 10000
a0_cutoff = 10
a = a0_cutoff*np.linspace(1,N_lattice_parameter,N_lattice_parameter)/N_lattice_parameter
e_embed = [rhofxn(v,[parameters]) for v in a]

import matplotlib.pyplot as plt

plt.plot(a,e_embed)
plt.show()
