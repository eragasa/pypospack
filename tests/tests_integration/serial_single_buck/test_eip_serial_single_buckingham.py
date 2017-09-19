import time
import pypospack.pyposmat as pyposmat
import pypospack.eipfit as eipfit



class TestEipSerialSingleBuckingham(object):
    def _set_up(self):
        # define the set of parameters
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

    def test_init(self):
        self.evaluator = eipfit.EipFittingEngine()
        

if __name__ == "__main__":
   
