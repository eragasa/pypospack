import pypospack.potential as potential
import numpy as np
from collections import OrderedDict

if __name__ == "__main__":
    symbols = ['Ni']
    param_names = ['NiNi_D0', 'NiNi_a', 'NiNi_r0']
    param_dict = OrderedDict()
    param_dict['NiNi_D0'] = 0.001114
    param_dict['NiNi_a'] = 3.429506
    param_dict['NiNi_r0'] = 2.6813
    morse = potential.MorsePotential(symbols=symbols)
    assert type(morse.potential_type) is str
    assert morse.potential_type == 'morse' 
    assert type(morse.symbols) is list
    assert morse.symbols == symbols
    assert type(morse.param_names) is list
    assert morse.param_names == param_names
    assert type(morse.is_charge) is bool

    if True:
        print("potential_type:{}".format(morse.potential_type))
        print("symbols:{}".format(morse.symbols))
        print("param_names:{}".format(morse.param_names))
        print("is_charge:{}".format(morse.is_charge))

    r = r_max * np.linespace(1,100,N_r)/100
    morse.evaluate(r,param_dict)

    symbols = ['Ni','Al']
    param_names = ['NiNi_D0', 'NiNi_a', 'NiNi_r0']
    morse = potential.MorsePotential(symbols=symbols)
    assert type(morse.potential_type) is str
    assert morse.potential_type == 'morse' 
    
    print("potential_type:{}".format(morse.potential_type))
    print("symbols:{}".format(morse.symbols))
    print("param_names:{}".format(morse.param_names))

if False:
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

            pair_plot = PairPotentialPlot(r,V) 
