import pytest
from collections import OrderedDict
import numpy as np

def create_eam_r_array(r_max,N_r):
    _r = r_max*np.linspace(1,N_r,N_r)/N_r
    return _r

def create_eam_rho_array(rho_max,N_rho):
    _rho = rho_max*np.linspace(1,N_rho,N_rho)/N/rho
    return _rho

def calculate_equilibrium_volume_per_atom(lattice_info):
    for k,v in lattice_info.items():
        if v['lattice_type'] == 'fcc':
            _a0 = v['lattice_parameter']
            _V = _a0**3

        v['volume_per_atom'] = _V

def calculate_equilibrium_interatomic_distance(lattice_info):
    for k,v in lattice_info.items():
        if v['lattice_type'] == 'fcc':
            _a0 = v['lattice_parameter']
            _re = 1/(2**0.5)*_a0

        v['equilbrium_interatomic_distance'] = _re

def convert_gpa_to_atom_units(pressure_in_gpa):
    _kg_to_g_per_mole = 1. #??
    _m_to_angs = 1.e10
    _s_to_ps = 1.e12
    # GPa in SI: N m^-2 == kg m s^-2 m^-2 = kg s^-2 m^-1
    P = pressure_in_gpa * 1.e9 # kg/m/s/s
    P = P * _kg_to_g_per_mole

from pypospack.potentials.eam_eos_foiles import func_eam_embed_foiles

from pypospack.potential import EamPotential

from pypospack.potential import BornMayerPotential
from pypospack.potential import ExponentialDensityFunction
from pypospack.potential import RoseEquationOfStateEmbeddingFunction

symbols = ['Ni']

func_pair_name = 'bornmayer'
func_density_name = 'eam_dens_exp'
func_embedding_name = 'eam_embed_eos_rose'

lattice_info = OrderedDict()
lattice_info['Ni'] = OrderedDict()
lattice_info['Ni']['equlibrium_interatomic_distance']
lattice_info['Ni']['cohesive_energy']
lattice_info['Ni']['bulk_modulus']
lattice_info['Ni']['lattice_parameter']
lattice_info['Ni']['bulk_modulus']

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

r_max = 10.
N_r = 1000

rho_max = 100000.
N_rho = 100
r = create_eam_r_array(r_max,N_r)
rho = create_eam_rho_array(rho_max,N_rho)

parameter_names_pair = ['p_NiNi_D0', 'p_NiNi_a', 'p_NiNi_r0']
parameter_names_density = ['d_Ni_rho0', 'd_Ni_beta', 'd_Ni_r0']
parameter_names_embedding = ['e_Ni_ecoh', 'e_Ni_latticetype', 'e_Ni_B', 'e_Ni_a0']
parameter_names = \
        parameter_names_pair\
        + parameter_names_density\
        + parameter_names_embedding

lattice_info = OrderedDict()
lattice_info['Ni'] = OrderedDict()
lattice_info['Ni']['lattice_type'] = 'fcc'
lattice_info['Ni']['cohesive_energy'] = -4.5
lattice_info['Ni']['bulk_modulus'] = 162
lattice_info['Ni']['lattice_parameter'] = 3.52


calculate_equilibrium_volume_per_atom(lattice_info)
calculate_equilibrium_interatomic_distance(lattice_info)

def test__init__():

    pot = EamPotential(symbols=symbols,
            func_pair=func_pair_name,
            func_density=func_density_name,
            func_embedding=func_embedding_name)

    assert pot.potential_type == 'eam'
    assert pot.is_charge is False

    assert pot.symbols == symbols

    assert type(pot.obj_pair) is BornMayerPotential
    assert type(pot.obj_embedding) is RoseEquationOfStateEmbeddingFunction
    assert type(pot.obj_density) is ExponentialDensityFunction

    assert pot.obj_pair.parameter_names == parameter_names_pair
    assert pot.obj_embedding.parameter_names == parameter_names_embedding
    assert pot.obj_density.parameter_names == parameter_names_density

def test__evaluate():

    pot = EamPotential(symbols=symbols,
            func_pair=func_pair_name,
            func_density=func_density_name,
            func_embedding=func_embedding_name)
    pot.evaluate(
            r = r,
            rho = rho,
            rcut = r_max,
            parameters = parameters)

