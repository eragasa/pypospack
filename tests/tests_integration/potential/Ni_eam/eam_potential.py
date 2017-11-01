import pyflamestk.potential as potential
import pyflamestk.eam as eam
import numpy as np

# configuration inforamtion
symbols = ['Ni']
name_eam_pair_func = 'morse'
name_eam_dens_func = 'exp'
name_eam_embed_func = 'universal'
name_eam_cutoff_func = ''
eam_pot = potential.EamPotential(symbols)

# setup functional forms
if name_eam_pair_func == 'morse':
    eam_pot.pair_function = potential.MorsePotential(symbols)

if name_eam_dens_func = 'exp':
    eam_pot.density_function = potential.ExponentialDensityFunction(symbols)

if name_eam_embed_func = 'universal':
    eam_pot.embedding_function = potential.UniversalEmbeddingFunction(symbols)

# cutoff type
if name_eam_cutoff_func == '':
    eam_pot.cutoff_function = ''
print(eam_pot.parameter_names)

# create parameter diction
param_dict = {}
param_dict['p.NiNi_D0'] = 1
param_dict['p.NiNi_a'] = 1
param_dict['p.NiNi_r0'] = 1
param_dict['d.Ni_rho0'] = 1
param_dict['d.Ni_beta'] = 1
param_dict['d.Ni_r0'] = 1
param_dict['e.Ni_F0'] = 1
param_dict['e.Ni_p'] = 1
param_dict['e.Ni_q'] = 1
param_dict['e.Ni_F1'] = 1
param_dict['rcut_g'] = 1
param_dict['rhocut_g'] = 1
param_dict['p.NiNi_hcut'] = 1
param_dict['d.Ni_hcut'] = 1

# additional
Nr = 500
Nrho = 500

rcut_g = param_dict['rcut_g']
rhocut_g = param_dict['rhocut_g']

# generate evenly spaced Nr points between r_low, r_high
# units in angstroms
r_low = 0
r_high = rcut_g
r = np.linspace(r_low,r_high,Nr)
dr = r[1]-r[0]

# generate evently spaced NRho points between rho_low, and rho_high
# units for rho is unclear, but should be a measure of electron density
rho_low = 0
rho_high = rhocut_g
rho = np.linspace(rho_low,rho_high,Nrho)
drho = rho[1]-rho[0]

dens_p_dict = {}
embed_p_dict = {}
for s in symbols:
    dens_p_dict[s] = {}
    embed_p_dict[s] = {}
    for k,v in param_dict.items():
        if k.startswith('e.{}'.format(s)):
            key = k.replace('e.','')
            embed_p_dict[s][key] = v
        elif k.startswith('d.{}'.format(s)):
            key = k.replace('d.','')
            dens_p_dict[s][key] = v

pair_p_dict = {}
N_symbols = len(symbols)
for i in range(N_symbols):
    for j in range(N_symbols):
        if i >= j:
            pair = '{}{}'.format(symbols[i],symbols[j])
            pair_p_dict[pair] = {}
            for k,v in param_dict.items():
                if k.startswith('p.{}'.format(pair)):
                    key = k.replace('p.','')
                    pair_p_dict[pair][key] = v

print(dens_p_dict)
print(embed_p_dict)
print(pair_p_dict)

