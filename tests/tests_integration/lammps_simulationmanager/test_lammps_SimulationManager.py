import copy
import pypospack.potfit as potfit
import pypospack.qoi as qoi
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

if __name__ == "__main__":
    qoi_config_file = "pypospack.qoi.yaml"
    structure_config_file = "pypospack.structure.yaml"
    potential_config_file = "pypospack.buckingham.yaml"

    qoi_info = potfit.QoiDatabase()
    qoi_info.read(qoi_config_file)

    structure_info = potfit.StructureDatabase()
    structure_info.read(structure_config_file)
    
    potential_info = potfit.PotentialInformation()
    potential_info.read(potential_config_file)

    qoi_manager = qoi.QoiManager(qoi_info)
    required_simulations = qoi_manager.get_required_simulations()
    for k,v in qoi_manager.required_simulations.items():
        print(k,'\n\t',v)
    sim_manager = lammps.SimulationManager()

    
