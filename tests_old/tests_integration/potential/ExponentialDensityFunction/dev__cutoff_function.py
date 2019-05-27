from collections import OrderedDict
from pypospack.potential import ExponentialDensityFunction
import matplotlib.pyplot as plt

r_max = 10
N_r = 100

import copy
import numpy as np
r = r_max * np.linspace(1,N_r,N_r)/N_r
if __name__ == "__main__":
    
    symbols = ['Ni']
    parameters = OrderedDict()
    parameters['Ni_rho0'] = 1 
    parameters['Ni_beta'] = 10
    parameters['Ni_r0'] = 10
    dens = ExponentialDensityFunction(symbols=symbols)
    
    dens_nocut = copy.deepcopy(dens.evaluate(r,parameters,r_cut=None))
    dens_rcut = copy.deepcopy(dens.evaluate(r,parameters,r_cut=4.))

    fig,ax = plt.subplots(1,1)
    ax.plot(r,dens_nocut['Ni'],color='b',label='nocut')
    ax.plot(r,dens_rcut['Ni'],color='k',label='cut')

    fig.savefig('Ni_exp_rcut.png')
