import pytest

import pypospack.potential as potential

class TestPotential(object):
    def setup__init_one_symbol(self):
        self.symbols = ['Ni']
        self.potential = potential.Potential(symbols=self.symbols)

    def test__init_one_symbol_w_list(self):
        self.symbols = ['Ni']
        self.potential = potential.Potential(symbols=self.symbols)

    def test__init_one_symbol_w_str(self):
        self.symbols = ['Ni']
        str_symbol = 'Ni'
        self.potential = potential.Potential(symbols='Ni')

    def test__init_one_symbol__attribute_symbol(self):
        self.setup__init_one_symbol()

        assert isinstance(self.potential.symbols, list)
        
        for i,v in enumerate(self.symbols):
            assert self.potential.symbols[i] == v

    def setup__init_two_symbols(self):
        self.symbols = ['Mg,O']
        self.potential = potential.Potential(symbols=self.symbols)

    def test__init_two_symbols__attribute_symbol(self):
        self.setup__init_two_symbols()

        assert isinstance(self.potential.symbols, list)

        for i,v, in enumerate(self.symbols):
            assert self.potential.symbols[i] == v

if __name__ == "__main__":

