import unittest
import pypospack.potential as potential

class test_potential_init_one_symbol(unittest.TestCase):

    def setUp(self):
        self.symbols = ['Ni']
        self.pot = potential.Potential(symbols = self.symbols)

    def test_attr_symbols(self):
        self.assertSequenceEqual(self.pot.symbols, self.symbols)

class test_potential_init_two_symbol(unittest.TestCase):
    def setUp(self):
        self.symbols = ['Mg','O']
        self.pot = potential.Potential(symbols = self.symbols)

    def test_attr_symbols(self):
        self.assertSequenceEqual(self.pot.symbols, self.symbols)

class test_buckingham_init(unittest.TestCase):
    def setUp(self):
        self.symbols = ['Mg','O']
        self.param_names = ['chrg_Mg', 'chrg_O', 
                            'MgMg_A', 'MgMg_rho', 'MgMg_C', 
                            'OO_A', 'OO_rho', 'OO_C', 
                            'MgO_A', 'MgO_rho', 'MgO_C']
        self.pot = potential.Buckingham(symbols = self.symbols)
    
    def test_attr_symbols(self):
        self.assertSequenceEqual(self.pot.symbols, self.symbols)

    def test_attr_param_names(self):
        self.assertSequenceEqual(self.pot.param_names, self.param_names)

    def test_attr_param_dict(self):
        keys = sorted(list(self.pot.param.keys()))
        param_names = sorted(self.param_names)
        self.assertSequenceEqual(keys,param_names)

class test_buckingham_clone(unittest.TestCase):
    def setUp(self):
        self.symbols = ['Mg','O']
        self.param_names = ['chrg_Mg', 'chrg_O', 
                            'MgMg_A', 'MgMg_rho', 'MgMg_C', 
                            'OO_A', 'OO_rho', 'OO_C', 
                            'MgO_A', 'MgO_rho', 'MgO_C']
        self.pot_orig = potential.Buckingham(symbols = self.symbols)
        self.pot_clone = self.pot_orig.copy()

    def test_attr_symbols(self):
        self.assertSequenceEqual(\
                self.pot_clone.symbols,
                self.pot_orig.symbols)

    def test_attr_param_names(self):
        self.assertSequenceEqual(\
                self.pot_clone.param_names,
                self.pot_orig.param_names)

    def test_attr_param_dict(self):
        for p in self.param_names:
            self.assertEqual(
                    self.pot_clone.param[p],
                    self.pot_orig.param[p])

if __name__ == '__main__':
    unittest.main()
    #symbols = ['Mg','O']
    #pot = potential.Buckingham(symbols=symbols)
    #print(pot.param_names)
    #print(pot.param)
