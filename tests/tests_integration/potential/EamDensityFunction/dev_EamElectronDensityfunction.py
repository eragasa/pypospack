import pypospack.potential as potential
from collections import OrderedDict
import numpy as np

import copy
import matplotlib.pyplot as plt
class EamElectronicDensityFunctionPlot(object):
    def __init__(self,r,dens,fname_plot='electronic_density.png'):
        assert isinstance(r,np.ndarray)
        assert isinstance(dens,np.ndarray)

        self.r = copy.deepcopy(r)
        self.dens = copy.deepcopy(dens)
        self.fname_plot=fname_plot

    def plot(self):
        self.figure, self.ax = plt.subplots(nrows=1,ncols=1)
        self.ax.plot(self.r,self.dens)

        self.figure.savefig(self.fname_plot)
        plt.close(self.figure)

symbols = ['Ni']
pot = potential.EamElectronDensityFunction(symbols=['Ni'])

attribute_dict = OrderedDict()
attribute_dict['symbols'] = pot.symbols
attribute_dict['potential_type'] = pot.potential_type
attribute_dict['param_names'] = pot.param_names

for k,v in attribute_dict.items():
    line = "{} = {}".format(k,v)
    print(line)

print(80*'-')
print('EXPONENTIAL DENSITY FUNCTION')
print(80*'-')

symbols = ['Ni']
pot = potential.ExponentialDensityFunction(symbols=symbols)

attribute_dict = OrderedDict()
attribute_dict['symbols'] = pot.symbols
attribute_dict['potential_type'] = pot.potential_type
attribute_dict['param_names'] = pot.param_names

for k,v in attribute_dict.items():
    line = "{} = {}".format(k,v)
    print(line)

r_min = 0
r_max = 10
N_r   = 100
r = np.linspace(r_min,r_max,N_r)

param_dict = {}
param_dict['d.Ni.rho0'] = 1
param_dict['d.Ni.beta'] = 1
param_dict['d.Ni.r0'] = 3
dens = pot.evaluate(r=r,symbol=symbols,params=param_dict,rcut=False)

dens_plot = EamElectronicDensityFunctionPlot(r,dens)
dens_plot.plot()

