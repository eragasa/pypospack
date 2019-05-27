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
Ni_eam_potential_definition['r_cut'] = 10.0
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

def test_1sym__write____morse_exponential_universal():
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
