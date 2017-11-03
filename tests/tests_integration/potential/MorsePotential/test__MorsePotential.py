import pytest
from collections import OrderedDict
import numpy as np
import pypospack.potential as potential

class Test__MorsePotential__one_element(object):

    def setup_variables(self):
        self.symbols = ['Ni']
        self.param_dict = OrderedDict()
        self.param_dict['NiNi_D0'] = 0.001114
        self.param_dict['NiNi_a'] = 3.429506
        self.param_dict['NiNi_r0'] = 2.6813

        self.r_max = 11.
        self.N_r = 500

    def test__init(self):
        self.setup_variables()
        self.morse = potential.MorsePotential(symbols=self.symbols)

    def test__evaluate(self):
        self.setup_variables()
        self.morse = potential.MorsePotential(symbols=self.symbols)
        r = self.r_max * np.linspace(1,100,self.N_r)/100
        self.morse.evaluate(r,self.param_dict)

        assert isinstance(self.morse.potential,OrderedDict)
        for pair_key,pot in self.morse.potential.items():
            assert isinstance(pot,np.ndarray)

class Test__MorsePotential__two_elements(object):

    def setup_variables(self):
        self.symbols = ['Ni','Al']
        self.symbol_pairs = [['Ni','Ni'],['Ni','Al'],['Al','Al']]
        self.param_dict = OrderedDict()
        self.param_dict['NiNi.D0'] = 0.001114
        self.param_dict['NiNi.a'] = 3.429506
        self.param_dict['NiNi.r0'] = 2.6813
        self.param_dict['NiAl.D0'] = 0.001114
        self.param_dict['NiAl.a'] = 3.429506
        self.param_dict['NiAl.r0'] = 2.6813
        self.param_dict['AlAl.D0'] = 0.001114
        self.param_dict['AlAl.a'] = 3.429506
        self.param_dict['AlAl.r0'] = 2.6813

    def test__init(self):
        self.setup_variables()
        
