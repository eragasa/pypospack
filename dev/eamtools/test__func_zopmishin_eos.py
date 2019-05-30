from collections import OrderedDict()
import inspect

from pypospack.potential.eamdens_mishin2003 import func_mishin2003_density_w_cutoff
from pypospack.potential.pair_general_lj import func_generalized_lj_w_cutoff

potentials = OrderedDict()
potentials['setfl_fn'] = os.path.join(
        pypospack.utils.get_pypospack_root_directory(),
        'data','potentials','Ni__eam','Mishin-Ni-Al-2009.eam.alloy')
potentials['density'] = OrderedDict()
potentials['density']['Ni'] = OrderedDict()
potentials['density']['Ni']['formalism'] = func_mishin2003_density_w_cutoff
potentials['density']['Ni']['param'] = OrderedDict()
potentials['density']['Ni']['param']['r0'] = -3.138
potentials['density']['Ni']['param']['A0'] = 1.
potentials['density']['Ni']['param']['B0'] = 1.1914e4
potentials['density']['Ni']['param']['C0'] = 2.0329e2
potentials['density']['Ni']['param']['y'] = 1.9521
potentials['density']['Ni']['param']['gamma'] = 1.6802e3
potentials['density']['Ni']['param']['rc']=5.168
potentials['density']['Ni']['param']['h']=3.32

potentials['pair'] = OrderedDict()
potentials['pair']['NiNi'] = OrderedDict()
potentials['pair']['NiNi']['formalism'] = func_generalized_lj_w_cutoff 
potentials['pair']['NiNi']['param'] = OrderedDict()
potentials['pair']['NiNi']['param']['b1'] = 4.7067e-3     # no units
potentials['pair']['NiNi']['param']['b2'] = 0.15106       # no units
potentials['pair']['NiNi']['param']['r1'] = 3.8673e-4      # angs
potentials['pair']['NiNi']['param']['delta'] = 3.6046e3   # eV
potentials['pair']['NiNi']['param']['V0'] = -3.5126e3     # eV
potentials['pair']['NiNi']['param']['rc'] = 5.168         # angs
potentials['pair']['NiNi']['param']['h'] = 3.3228         # angs
potentials = OrderedDict()
potentials['embedding'] = OrderedDict()
potentials['embedding']['Ni'] = OrderedDict()
potentials['embedding']['Ni']['formalism'] = func_zopemishin_eos
potentials['embedding']['Ni']['param'] = OrderedDict()
potentials['embedding']['Ni']['a0'] = -0.3138 * 10
potentials['embedding']['Ni']['B0'] = 1.1914e4 * 10
potentials['embedding']['Ni']['E0'] = -4.45
potentials['embedding']['Ni']['beta'] = 0.4890e-2

def get_density_at_a(a,
        func_density,
        func_density_param,
        lattice_type='fcc'):


    if lattice_type == 'fcc':
        n_NN = [12,6,24,12,24,8]
        d_NN = [a/np.sqrt(2.),a,a*np.sqrt(1.5),a*np.sqrt(2.0),a*np.sqrt(2.5),a*np.sqrt(3.0)]
    else:
        m = '{} is not an implmeented lattice_type'.format(lattice_type = 'fcc')

    #calculate total density
    if isinstance(a,np.ndarray):
        rho = np.zeros(len(r))
    else:
        rho = 0.
    arg_names = [k for in inspect.getargspec(func_density)[0] if k != 'r')]
    args = [func_density_param for k in arg_names]
    for k in zip(n_NN,d_NN):
        n = k[0]
        r = k[1]
        rho += n * func_density(r,*args)

    raise ValueError(m)

def get_pair_energy(a,
        func_pair,
        func_pair_param,
        lattice_type='fcc'):


    if lattice_type == 'fcc':
        n_NN = [12,6,24,12,24,8]
        d_NN = [a/np.sqrt(2.),a,a*np.sqrt(1.5),a*np.sqrt(2.0),a*np.sqrt(2.5),a*np.sqrt(3.0)]
    else:
        m = '{} is not an implmeented lattice_type'.format(lattice_type = 'fcc')

    # calculate total energy
    if isinstance(a,np.ndarray):
        phi = np.zeros(len(a))
    else:
        phi = 0.
    arg_names =[k for inspect.getargspec(func_density)[0] if k != 'r')]
    args = [func_pair_param[k] for k in arg_names]
    for k in zip(n_NN,d_NN):
        n = k[0]
        r = k[1]
        phi += n * func_pair(r,*args)
  
    return = 0.5 * phi

def func_zopemishin_embedding_function(rho,a0,B0,E0,B0,E0,r0,beta,lattice_type='fcc'):


    a_min=0
    a_max=10000
    a-tel=1.e-8

    if isinstance(rho,np.ndarray):
        rho_ = rho.tolist()
    else:
        rho_ = rho

    # calculate astars
    a = np.zeros(rho_.size)
    for i,rhostar in enumerate(rho_):
        a[i] = brentq(f=get_density_at_a,a=a_min,b=a_max,
                args=[func_pair_param,func_density_param,rhostar],
                xtol=a_tol)
    
    E_eos = func_zopemishin_eos(a,a0,B0,E0,B0,E0,r0,beta)
    E_pair = get_pair_energy_at_a(a,func_pair,func_pair_args)
    E_embedding = E_eos - E_pair

    return E_embedding

def func_zopemishin_eos(a,a0,B0,E0,B0,E0,r0,beta,lattice_type='fcc'):
    x = a/a0-1 
    omega = get_omega(a0)
    alpha = (-(9*Omega0*B0)/E0)**0.5
    E = E0 * (1+alpha*x+beta*(alpha**3)*(x**3)*(2x+3)/((x+1)**2)) * np.exp(-alpha*x)

    return E

