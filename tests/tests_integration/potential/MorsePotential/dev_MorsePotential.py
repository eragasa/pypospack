import pypospack.potential as potential
import numpy as np
from collections import OrderedDict

def print_potential_attributes(pot):
    attribute_dict = OrderedDict()
    attribute_dict['symbols'] = pot.symbols
    attribute_dict['potential_type'] = pot.potential_type
    attribute_dict['param_names'] = pot.param_names
    attribute_dict['is_charge'] = pot.is_charge 
    for k,v in attribute_dict.items():
        print("{} = {}".format(k,v))

import copy
import matplotlib.pyplot as plt

class PairPotentialPlot(object):
    def __init__(self,r,V,fname_plot='pair_potential.png'):
        assert isinstance(r,np.ndarray)
        assert isinstance(V,np.ndarray)

        self.r = copy.deepcopy(r)
        self.V = copy.deepcopy(V)
        self.fname_plot=fname_plot

        self.v_lim = [-0.01,0.01]
        self.r_lim = [min(self.r),max(self.r)]

    def plot(self):
        self.figure, self.ax = plt.subplots(nrows=1,ncols=1)

        self.ax.plot(self.r,self.V)
        self.ax.set_xlim(self.r_lim)
        self.ax.set_ylim(self.v_lim)

        self.figure.savefig(self.fname_plot)
        plt.close(self.figure)

print(80*'-')
print('MorsePotential, one_symbol')
print(80*'-')

symbols = ['Ni']
param_dict_NiNi = OrderedDict()
param_dict_NiNi['NiNi.D0'] = 0.001114 #eV 
param_dict_NiNi['NiNi.a'] = 3.429506 #Angs^(-2)
param_dict_NiNi['NiNi.r0'] = 2.6813 #Angs

pot = potential.MorsePotential(symbols=symbols)
print_potential_attributes(pot)

#test_evaluate
r_min = 0
r_max = 10
N_r = 500
d_r = (r_max-r_min)/(N_r-1)
r = np.linspace(r_min,r_max,N_r)

#print('r:{},{},{},{}'.format(r_min,r_max,N_r,d_r))
#print('N_r_theo:{}'.format(N_r))
#print('N_r_actual:{}'.format(len(r)))
#print('d_r_theo:{}'.format(d_r))
#print('d_r_actual:{}'.format(r[1]-r[0]))

V = pot.evaluate(r,['Ni','Ni'],param_dict_NiNi,rcut=False)

pair_plot = PairPotentialPlot(r,V) 
pair_plot.plot()


print
print(80*'-')
print('MorsePotential_two_symbols')
print(80*'-')
symbols = ['Ni','Al']

pot = potential.MorsePotential(symbols=symbols)
print_potential_attributes(pot)

