import numpy as np

# functions for the mishin2003 density function
from pypospack.potential.eamdens_mishin2003 import func_cutoff_mishin2003
from pypospack.potential.eamdens_mishin2003 import func_mishin2003_density
from pypospack.potential.eamdens_mishin2003 import func_mishin2003_density_w_cutoff

def create_r(rmax,n):
    r = rmax/n * np.linspace(1,n,n)
    return r

from collections import OrderedDict

def determine_eam_pairs(symbols):
    """ determines the order of eam pairs

    the order for eam pairs is different than for the rest of pypospack

    Args:
        symbols (list): a list of symbols

    Returns:
        (list): a list of tuples
    """


    for i1, s1 in enumerate(symbols):
        for i2, s2 in enumerate(symbols):
            if i1 >= i2:


class EamPotentialFitter(object):

    def __init__(symbols):


Ni_dens_params = OrderedDict()
Ni_dens_params['r0'] = -3.138
Ni_dens_params['A0'] = 1.
Ni_dens_params['B0'] = 1.1914e4
Ni_dens_params['C0'] = 2.0329e2
Ni_dens_params['y'] = 1.9521
Ni_dens_params['gamma'] = 1.6802e3
Ni_dens_params['rc']=5.168
Ni_dens_params['h']=3.32

import matplotlib.pyplot as plt
def plot_fitted_density_functions(r,rho_file,rho_fit):
    fig,ax = plt.subplots(1,1)

    ax.plot(r,rho_file,label='from file')
    ax.plot(r,rho_fit,label='from_fit')
    ax.legend()
    plt.show()

def fit_density_functions(setfl_fn):
    from pypospack.eamtools import SeatonSetflReader
    setfl = SeatonSetflReader(path=setfl_fn)
    setfl.read()

    rdata = create_r(setfl.n_r*setfl.d_r,setfl.n_r)
    rhodata = setfl.density_function('Ni')
   
    param_names = ['r0','A0','B0','C0','y','gamma','rc','h']
    p0 = [Ni_dens_params[k] for k in param_names]

    from scipy.optimize import curve_fit
    popt,pcov = curve_fit(
            func_mishin2003_density_w_cutoff,
            rdata,rhodata,
            method='dogbox',
            p0=p0)
   
    plot_fitted_density_functions(
            rdata,
            rhodata,
            func_mishin2003_density_w_cutoff(rdata,*popt))
    for k in zip(param_names,p0,popt):
        print(k)

    return OrderedDict([k for k in zip(param_names,popt)])

from pypospack.potential.pair_general_lj import function_generalized_lj_pair

NiNi_pair_params = OrderedDict()
NiNi_pair_params['b1'] = 4.7067e-3     # no units
NiNi_pair_params['b2'] = 0.15106       # no units
NiNi_pair_params['r1'] = 3.8673e-4      # angs
NiNi_pair_params['delta'] = 3.6046e3   # eV
NiNi_pair_params['V0'] = -3.5126e3     # eV
NiNi_pair_params['rc'] = 5.168         # angs
NiNi_pair_params['h'] = 3.3228         # angs

def get_1NN_fcc(a0=3.52):
    return a0*(2**.5)/2

def func_generalized_lj_w_cutoff(r,b1,b2,r1,V0,delta,rc,h):

    phi = function_generalized_lj_pair(r,b1,b2,r1,V0,delta)
    psi = func_cutoff_mishin2003(r,rc,h)

    return psi*phi

def plot_pair_potential(r,phi_file,phi_fit):
    fig, ax = plt.subplots(1,1)

    ax.plot(r,phi_file,label='file')
    ax.plot(r,phi_fit,label='fit')
    ax.legend()
    ax.set_xlim(1.5,4)
    ax.set_ylim(-1,5)
    plt.show()
def fit_pair_potentials(setfl_fn):
    from pypospack.eamtools import SeatonSetflReader
    setfl = SeatonSetflReader(path=setfl_fn)
    setfl.read()


    param_names = ['b1','b2','r1','V0','delta','rc','h']
    param_0 = [NiNi_pair_params[k] for k in param_names]

    rdata = create_r(setfl.n_r*setfl.d_r,setfl.n_r)
    phidata = np.array(setfl.pair_function('NiNi'))/rdata
    rdata_2 = rdata[rdata>1.5]
    phidata_2 = phidata[rdata>1.5]

    param_0
    from scipy.optimize import curve_fit
    while True:
        param_opt,pocv = curve_fit(
                func_generalized_lj_w_cutoff,
                rdata_2,phidata_2,
                method='trf',
                p0=param_0)
       
        for k in zip(param_names,param_0,param_opt):
            print(k)

        if all([np.abs(k[1]/k[0]-1) < 0.01 for k in zip(param_opt,param_0)]):
            break
        param_0 = param_opt
    plot_pair_potential(
            rdata,
            phidata,
            func_generalized_lj_w_cutoff(rdata,*param_opt))

if __name__ == "__main__":
    import os
    import pypospack.utils
    from pypospack.potential import EamPotential

    from collections import OrderedDict

    setfl_filename = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'data','potentials','Ni__eam','Mishin-Ni-Al-2009.eam.alloy')

    print(get_1NN_fcc())
    fit_density_functions(setfl_fn=setfl_filename)
    fit_pair_potentials(setfl_fn=setfl_filename)
    exit()
