import unittest
import pypospack.potential as potential

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


class test_buckingham_to_string(unittest.TestCase):
    def setUp(self):
        self.symbols = ['Mg','O']
        self.param_names = ['chrg_Mg', 'chrg_O', 
                            'MgMg_A', 'MgMg_rho', 'MgMg_C', 
                            'OO_A', 'OO_rho', 'OO_C', 
                            'MgO_A', 'MgO_rho', 'MgO_C']

        self.param_dict = {}
        self.param_dict['chrg_Mg'] = +2.0
        self.param_dict['chrg_O']  = -2.0
        self.param_dict['MgMg_A']   = 0.0 
        self.param_dict['MgMg_rho'] = 0.5
        self.param_dict['MgMg_C']   = 0.0
        self.param_dict['MgO_A']    = 821.6
        self.param_dict['MgO_rho']  = 0.3242
        self.param_dict['MgO_C']    = 0.0
        self.param_dict['OO_A']     = 2274.00 
        self.param_dict['OO_rho']   = 0.1490
        self.param_dict['OO_C']     = 27.88
        self.pot_orig = potential.Buckingham(symbols = self.symbols)

    def test_to_string_w_param_dict(self):
        param_dict = {}
        param_dict['chrg_Mg'] = +2.0
        param_dict['chrg_O']  = -2.0
        param_dict['MgMg_A']   = 0.0 
        param_dict['MgMg_rho'] = 0.5
        param_dict['MgMg_C']   = 0.0
        param_dict['MgO_A']    = 821.6
        param_dict['MgO_rho']  = 0.3242
        param_dict['MgO_C']    = 0.0
        param_dict['OO_A']     = 2274.00 
        param_dict['OO_rho']   = 0.1490
        param_dict['OO_C']     = 27.88

        pot_orig = potential.Buckingham(symbols = symbols)
        print(pot_orig.to_string(param_dict))
        

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
    # unittest.main()

    
    symbols = ['Mg','O']
    param_names = ['chrg_Mg', 'chrg_O', 
                        'MgMg_A', 'MgMg_rho', 'MgMg_C', 
                        'OO_A', 'OO_rho', 'OO_C', 
                        'MgO_A', 'MgO_rho', 'MgO_C']

    param_dict = {}
    param_dict['chrg_Mg'] = +2.0
    param_dict['chrg_O']  = -2.0
    param_dict['MgMg_A']   = 0.0 
    param_dict['MgMg_rho'] = 0.5
    param_dict['MgMg_C']   = 0.0
    param_dict['MgO_A']    = 821.6
    param_dict['MgO_rho']  = 0.3242
    param_dict['MgO_C']    = 0.0
    param_dict['OO_A']     = 2274.00 
    param_dict['OO_rho']   = 0.1490
    param_dict['OO_C']     = 27.88
    pot_orig = potential.Buckingham(symbols = symbols)
    print(pot_orig.to_string(param_dict))
