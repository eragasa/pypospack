import pytest
import os
from collections import OrderedDict

# <----- TESTING HARNESS VARIABLES

Ni_eam_potential_definition = OrderedDict()
Ni_eam_potential_definition['potential_type'] = 'eam'
Ni_eam_potential_definition['setfl_filename']=None
Ni_eam_potential_definition['pair_type']='morse'
Ni_eam_potential_definition['density_type']='eam_dens_exp'
Ni_eam_potential_definition['embedding_type']='eam_embed_universal'
Ni_eam_potential_definition['N_r'] = 2000
Ni_eam_potential_definition['r_max'] = 10.0
Ni_eam_potential_definition['r_cut'] = 8.9
Ni_eam_potential_definition['N_rho'] = 2000
Ni_eam_potential_definition['rho_max'] = 10.0
Ni_eam_potential_definition['symbols'] = ['Ni']

Ni_eam_parameters = OrderedDict()
Ni_eam_parameters['p_NiNi_D0'] = 0.001114
Ni_eam_parameters['p_NiNi_a'] = 3.429506
Ni_eam_parameters['p_NiNi_r0'] = 2.6813
Ni_eam_parameters['d_Ni_rho0'] = 10.0
Ni_eam_parameters['d_Ni_beta'] = 5.0
Ni_eam_parameters['d_Ni_r0'] = 2.0
Ni_eam_parameters['e_Ni_F0'] = 4.10341782e-3
Ni_eam_parameters['e_Ni_p'] = 8.96274624
Ni_eam_parameters['e_Ni_q'] = 8.95940869
Ni_eam_parameters['e_Ni_F1'] = -3.09

configuration = OrderedDict()
configuration['potential'] = Ni_eam_potential_definition
configuration['parameters'] = Ni_eam_parameters

#<------------- unpack dictionary
symbols = configuration['potential']['symbols']
func_pair = configuration['potential']['pair_type']
func_density = configuration['potential']['density_type']
func_embedding = configuration['potential']['embedding_type']

parameters = configuration['parameters']
#<------------ setup for testing
from pypospack.potential import EamPotential
eam = EamPotential(
        symbols=symbols,
        func_pair=func_pair,
        func_density=func_density,
        func_embedding=func_embedding)

a0 = 3.50803
sites = ['T','O','1NN','2NN','3NN']
N = OrderedDict()
N['T'] = 4
N['O'] = 4
N['1NN'] = 12
N['2NN'] = 6
N['3NN'] = 24
da = OrderedDict()
da['T'] = 0.433
da['O'] = 0.866
da['1NN']= 0.707
da['2NN']= 1.00
da['3NN'] = 1.225

rcut = 0.5 * (da['2NN'] + da['3NN'])
rho = OrderedDict()
for s in sites:
    _pair = eam.evaluate_density(
            r=a0*da[s],
            parameters=parameters)
    for p in _pair:
        if p not in rho: rho[p] = OrderedDict()
        rho[p][s] = _pair[p]

print(rho)
_pair = [p for p in rho]
for p in _pair:
    for s in sites:
        print("{s:^10}{N_s:^10}{da:^10.4f}{rho:^15.4e}{ttl_rho:^15.4e}".format(
            s=s,
            N_s=N[s],
            da=da[s],
            rho=rho[p][s],
            ttl_rho=N[s]*rho[p][s]))
sites_lte_1NN = ['1NN']
rho_lte_1NN = 0.


filename = "Ni.eam.alloy"
Nr = configuration['potential']['N_r']
rmax = configuration['potential']['r_max']
rcut = configuration['potential']['r_cut']
Nrho = configuration['potential']['N_rho']
rhomax = configuration['potential']['rho_max']
eam.write_setfl_file(
        filename=filename,
        symbols=symbols,
        Nr=Nr,
        rmax=rmax,
        rcut=rcut,
        Nrho=Nrho,
        rhomax=rhomax,
        parameters=parameters)

assert os.path.isfile(filename)
