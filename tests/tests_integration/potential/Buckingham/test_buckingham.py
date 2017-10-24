import pytest

class TestBuckingham_init(object):

    def setup_init(self):
        self.symbols = ['Mg','O']
        self.is_charge = True
        self.param_names = [
            'chrg_Mg','chrg_O',
            'MgMg_A','MgMg_C','MgMg_rho',
            'MgO_A','MgO_C','MgO_rho',
            'OO_A','OO_C','OO_rho']

        import pypospack.potential as potential
        self.potential = potential.Buckingham(symbols=self.symbols)

        assert isinstance(self.potential,potential.Buckingham)

    # test that inherited attributes are established correctly
    def test_init__attribute__potential_type(self):
        self.setup_init()
        assert self.potential.potential_type == 'buckingham'
        
    def test_init__attribute__symbols(self):
        self.setup_init()

        for i,v in enumerate(self.potential.symbols):
            assert v == self.symbols[i]

    def test_init__attribute__param_names(self):
        self.setup_init()
        
        for v in self.param_names:
            assert v in self.potential.param_names

    def test_init__attribute__is_charge(self):
        self.setup_init()

        assert self.potential.is_charge is self.is_charge

    def test_init__attribute__param(self):
        self.setup_init()

        assert isinstance(self.potential.param,dict)
