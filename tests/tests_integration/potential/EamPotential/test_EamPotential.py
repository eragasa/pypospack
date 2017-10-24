import pytest

class Test_EamPotential_one_symbol(object):
    def setup_eam_potential(self):
        import pypospack.potential as potential
        self.symbols = ['Ni']
        self.pot = potential.EamPotential(symbols=self.symbols)

    def test__init(self):
        self.setup_eam_potential()

    def test__init__attribute_symbols(self):
        self.setup_eam_potential()

        assert isinstance(self.pot.symbols,list)
        assert len(self.pot.symbols) == len(self.symbols)

        for i,v, in enumerate(self.symbols):
            assert self.pot.symbols[i] == v

class Test_EamPotential_two_symbol(object):
    def setup_eam_potential(self):
        import pypospack.potential as potential
        self.symbols = ['Ni','Al']
        self.pot = potential.EamPotential(symbols=self.symbols)

    def test__init(self):
        self.setup_eam_potential()

    def test__init__attribute_symbols(self):
        self.setup_eam_potential()

        assert isinstance(self.pot.symbols,list)
        assert len(self.pot.symbols) == len(self.symbols)

        for i,v, in enumerate(self.symbols):
            assert self.pot.symbols[i] == v
