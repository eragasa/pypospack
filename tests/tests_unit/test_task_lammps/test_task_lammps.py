import unittest
import pypospack.potential as potential
import pypospack.task.lammps as lammps_task

class test_task_LammpsMinimizeStructure_buckingham(unittest.TestCase):

    def setUp(self):
        symbols = ['Mg','O']
        fname = 'potential.mod'

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

        self.buck = potential.Buckingham(symbols = symbols)
        self.buck.param = param_dict.copy()
        self.task = lammps_task.LammpsMinimizeStructure(fname, self.buck)

    
    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()


