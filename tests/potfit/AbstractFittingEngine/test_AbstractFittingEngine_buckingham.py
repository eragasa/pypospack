import copy
import pypospack.potfit as potfit
import pypospack.potential as potential
import pypospack.lammps as lammps
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

    print(80*'-')
    print('\t Inspecting SimulationManager')
    sim_mgr = fitter.simulation_manager
    sim_mgr.variable_dict= {}
    for task_name, obj_lammps_task in sim_mgr.obj_lammps_tasks.items():
        for var_name, var_value in obj_lammps_task['obj'].results.items():
            k = '{}.{}'.format(task_name,var_name)
            sim_mgr.variable_dict[k] = var_value

    assert isinstance(sim_mgr, lammps.SimulationManager)
    # assert isinstance(sim_mgr.variable_names, list)
    assert isinstance(sim_mgr.variable_dict, dict)
    for k,v in sim_mgr.variable_dict.items():
        print(k,v)

    print(80*'-')
    print('\t Inspecting the QoiManager')
    for qoi_name,qoi in fitter.qoi_manager.obj_qois.items():
        print(qoi_name,qoi)
        qoi.calculate_qoi(fitter.variables)
        

    print(80*'-')
    print('\t Inspecting the QoiDatabase')
    for k,v, in fitter.qoi_info.qois.items():
        print(k,v)

    print(80*'-')
    print('\t Inspecting the SimulationManager')
    for k,v in fitter.simulation_manager.obj_lammps_tasks.items():
        print(k,type(v['obj']))
        #for var_name, var_value in v['variables'].items():
        #    print('\t',var_name,var_value)

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

   
