from collections import OrderedDict
from pypospack.potential import MorsePotential
import matplotlib.pyplot as plt

r_max = 10
N_r = 100

import copy
import numpy as np
r = r_max * np.linspace(1,N_r,N_r)/N_r
if __name__ == "__main__":
    
    symbols = ['Ni']
    parameters = OrderedDict()
    parameters['NiNi_D0'] = 0.001114
    parameters['NiNi_a'] = 3.429506
    parameters['NiNi_r0'] = 2.6813
    
    pair = MorsePotential(symbols=symbols)
    
    pair_nocut = copy.deepcopy(pair.evaluate(r,parameters,r_cut=None))
    pair_rcut = copy.deepcopy(pair.evaluate(r,parameters,r_cut=4.))

    fig,ax = plt.subplots(1,1)
    ax.plot(r,pair_nocut['NiNi'],color='b',label='nocut')
    ax.plot(r,pair_rcut['NiNi'],color='k',label='cut')

    ax.set_ylim([-.001,.001])
    fig.savefig('Ni_morse_rcut.png')
