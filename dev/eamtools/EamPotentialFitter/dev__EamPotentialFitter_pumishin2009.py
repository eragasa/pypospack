import os,inspect
from collections import OrderedDict

import numpy as np
from scipy.optimize import brentq

import pypospack.utils
from pypospack.eamtools import create_r

from pypospack.potential.eamdens_mishin2004 import func_density_mishin2004_w_cutoff
from pypospack.potential.pair_general_lj import func_pair_generalized_lj_w_cutoff
from pypospack.potential.eamembed_eos_zopemishin import func_zopemishin_eos
from pypospack.potential.eamembed_eos_zopemishin import func_zopemishin_embedding_function
# currently in development
from potentialfitter import EamPotentialFitter
from eossolver import get_pair_energy_at_a
from eossolver import get_density_at_a 

# Mishin2003 -- original parameterization
potentials = OrderedDict()
potentials['setfl_fn'] = os.path.join(
        pypospack.utils.get_pypospack_root_directory(),
        'data','potentials','Ni__eam','Mishin-Ni-Al-2009.eam.alloy')
potentials['density'] = OrderedDict()
potentials['density']['Ni'] = OrderedDict()
potentials['density']['Ni']['formalism'] = func_density_mishin2004_w_cutoff
potentials['density']['Ni']['param'] = OrderedDict()
potentials['density']['Ni']['param']['r0'] = -3.138
potentials['density']['Ni']['param']['A0'] = 1.
potentials['density']['Ni']['param']['B0'] = 1.1914e4
potentials['density']['Ni']['param']['C0'] = 2.0329e2
potentials['density']['Ni']['param']['y'] = 1.9521
potentials['density']['Ni']['param']['gamma'] = 1.6802e3
potentials['density']['Ni']['param']['rc']=5.168
potentials['density']['Ni']['param']['hc']=3.32
potentials['density']['Ni']['param']['h0']=0.75

potentials['pair'] = OrderedDict()
potentials['pair']['NiNi'] = OrderedDict()
potentials['pair']['NiNi']['formalism'] = func_pair_generalized_lj_w_cutoff 
potentials['pair']['NiNi']['pair'] = ['Ni', 'Ni']
potentials['pair']['NiNi']['param'] = OrderedDict()
potentials['pair']['NiNi']['param']['b1'] = 4.7067e-3     # no units
potentials['pair']['NiNi']['param']['b2'] = 0.15106       # no units
potentials['pair']['NiNi']['param']['r1'] = 3.8673e-4      # angs
potentials['pair']['NiNi']['param']['delta'] = 3.6046e3   # eV
potentials['pair']['NiNi']['param']['V0'] = -3.5126e3     # eV
potentials['pair']['NiNi']['param']['rc'] = 5.168         # angs
potentials['pair']['NiNi']['param']['hc'] = 3.3228         # angs
potentials['pair']['NiNi']['param']['h0'] = 3.3228         # angs

potentials['embedding'] = OrderedDict()
potentials['embedding']['Ni'] = OrderedDict()
potentials['embedding']['Ni']['formalism'] = func_zopemishin_eos
potentials['embedding']['Ni']['param'] = OrderedDict()
potentials['embedding']['Ni']['param']['a0'] = 0.3138 * 10
potentials['embedding']['Ni']['param']['B'] = 1.1914e4 * 10
potentials['embedding']['Ni']['param']['E0'] = -4.45
potentials['embedding']['Ni']['param']['beta'] = 0.4890e-2

setfl_fn = os.path.join(
        pypospack.utils.get_pypospack_root_directory(),
        'data','potentials','Ni__eam','Mishin-Ni-Al-2009.eam.alloy')
symbols = ['Ni']
o = EamPotentialFitter(symbols)
o.read_setfl_file(filename=setfl_fn)

print('fitting the pair potential')
o.fit_potential_pair(
        func_pair_potential=potentials['pair']['NiNi']['formalism'],
        symbol_pair=['Ni','Ni'],
        param0=potentials['pair']['NiNi']['param'],
        rlow=1.5)
assert o.formalisms['pair']['NiNi'] == potentials['pair']['NiNi']['formalism']
print(o.parameters['p0']['pair']['NiNi'])
print(o.parameters['popt']['pair']['NiNi'])

print('fitting the density function')
o.fit_density_function(
        func_density=potentials['density']['Ni']['formalism'],
        symbol='Ni',
        param0=potentials['density']['Ni']['param'],
        rlow=1.5)
assert o.formalisms['density']['Ni'] == potentials['density']['Ni']['formalism']
print(o.parameters['p0']['density']['Ni'])
print(o.parameters['popt']['density']['Ni'])

def plot_results(filename='fit_eam_punmishin.eps'):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(3,1)

    rmax = o.setfl_reader.n_r * o.setfl_reader.d_r
    rN = o.setfl_reader.n_r
    r = create_r(rmax,rN)
    ax[0].plot(
            r,
            o.setfl_reader.pair_function('NiNi'),
            label='from_file')
    ax[0].plot(
            r,
            o.formalisms['pair']['NiNi'](r,**o.parameters['popt']['pair']['NiNi']),
            label='fitted')

    ax[1].plot(
            r,
            o.setfl_reader.density_function('Ni'),
            label='from_file')
    ax[1].plot(
            r,
            o.formalisms['density']['Ni'](r,**o.parameters['popt']['density']['Ni']),
            label='fitted')
    
    ax[2].plot(
            rho,
            o.setfl_reader.embedding_function('Ni'),
            label='from_file')
    #ax[2].plot(
    #        rho,
    #        o.formalisms['embedding']['Ni'](r,**o.parameters['popt']['embedding']['Ni']),
    #        label='fitted')
    plt.show()
    fig.savefig(filename)

def func_zopemishin_embedding_function(rho,a0,B0,E0,beta,lattice_type='fcc'):
    """ fits the embedding function to the zope mishin equation of state

    This function only exists as a implementation prototype for optimization routines
    which requires an encapsulated function.  The EosFitter is a more general implementation,
    which should be used in most cases.

    Args:
        rho (numpy.ndarray): a list of densities for the embedding function to be evaluated at
        a0 (float): the equilibrium lattice parameters
        B0 (float): the equilibrium bulk modulus.
        beta (float): a shape parameter
        lattice_type (str): the type of lattice parameter, currently only 'fcc' is implemented

    Returns:
        (numpy.ndarray): evaluations of the embedding function for rho

    """

    func_density = o.formalisms['density']['Ni']
    func_density_param = o.parameters['popt']['density']['Ni']

    func_pair = o.formalisms['pair']['NiNi']
    func_pair_param = o.parameters['popt']['pair']['NiNi']

    func_eos = o.formalisms['embedding']['Ni']

    a_min=0
    a_max=10000
    a_tol=1.e-8

    if isinstance(rho,np.ndarray):
        rho_ = rho.tolist()
    else:
        rho_ = rho

    # calculate astars
    a = np.zeros(len(rho_))
    for i,rhostar in enumerate(rho_):
        a[i] = brentq(f=get_density_at_a,a=a_min,b=a_max,
                args=(rhostar,func_density,func_density_param),
                xtol=a_tol)
    
    E_eos = func_zopemishin_eos(a,a0,B0,E0,beta)
    E_pair = get_pair_energy_at_a(a,func_pair,func_pair_param)
    E_embedding = E_eos - E_pair

    return E_embedding

rhomax = o.setfl_reader.n_rho*o.setfl_reader.d_rho
rhoN = o.setfl_reader.n_rho
rho = create_r(rhomax,rhoN)

arg_names = [k for k in inspect.getargspec(func_zopemishin_embedding_function)[0] if k not in ['rho','lattice_type']]
print(arg_names)
args = [potentials['embedding']['Ni']['param'][k] for k in arg_names]

plot_results()
embedding = func_zopemishin_embedding_function(rho,*args)
assert embedding.size == rho.size
o.fit_eos_embedding_function(
        func_embedding=func_zopemishin_embedding_function,
        symbol='Ni',
        param0=potentials['embedding']['Ni']['param'],
        bounds=[

            (3., 1e4,     -6,      -0.4890e-1),
            (4,  1e6,     -4,       0.4890e-1)
            ]
        )
print(o.parameters['p0']['embedding']['Ni'])
print(o.parameters['popt']['embedding']['Ni'])
