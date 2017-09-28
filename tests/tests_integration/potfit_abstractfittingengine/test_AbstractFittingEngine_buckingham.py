import copy
import pypospack.potfit as potfit
import pypospack.potential as potential

def get_lewis_catlow_parameters():
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
    return copy.deepcopy(param_dict)

class TestAbstractFittingEngineBuckingham(object):

    def setup(self):
        self.potfit = potfit.AbstractFittingEngine(\
                fname_config_potential = 'pypospack.buckingham.yaml',
                fname_config_qoi = 'pypospack.qoi.yaml',
                fname_config_structures = 'pypospack.structure.yaml')

    def test_config(self):
        self.setup()
        assert isinstance(self.potfit, potfit.AbstractFittingEngine)
        assert isinstance(self.potfit.structure_info, potfit.StructureDatabase)
        assert isinstance(self.potfit.qoi_info, potfit.QoiDatabase)
        assert isinstance(self.potfit.potential_info, potfit.PotentialInformation)
        assert isinstance(self.potfit.obj_potential,potential.Buckingham)

    def test_config_qoi_names(self):
        self.setup()
        print(self.potfit.qoi_names)
        assert isinstance(self.potfit.qoi_names, list)

    def test_config_parameter_names(self):
        self.setup()
        print(self.potfit.parameter_names)
        assert isinstance(self.potfit.parameter_names,list)

    def test_config_free_parameters(self):
        self.setup()
        print(self.potfit.free_parameters)
        assert isinstance(self.potfit.free_parameters,list)

    def test_configure_potential(self):
        self.setup()
if __name__ == "__main__":
    fitter = potfit.AbstractFittingEngine(\
                fname_config_potential = 'pypospack.buckingham.yaml',
                fname_config_qoi = 'pypospack.qoi.yaml',
                fname_config_structures = 'pypospack.structure.yaml')

    fitter.evaluate_parameter_set(\
            param_dict = get_lewis_catlow_parameters())

    def print_qoi_names(fitter):
        assert isinstance(fitter, potfit.fitter)
        for n in fitter.qoi_names:
            print(n)

    def print_parameter_names(fitter):
        assert isinstance(fitter, potfit.AbstractFittingEngine)
        for n in fitter.parameter_names:
            print(n)

    def print_free_parameter_names(fitter):
        assert isinstance(potfit.AbstractFittingEngine,fitter)
        for n in fitter.free_parameters:
            print(p)

    # TODO: move this the correct interation test
    def print_free_parameter_names__potential_info(pot_info):
        assert isinstance(potfit.PotentialInfo,pot_info)
        for n in pot_info.free_parameters:
            print(p)

   
