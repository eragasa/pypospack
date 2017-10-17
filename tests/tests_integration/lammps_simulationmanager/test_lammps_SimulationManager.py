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
    is_debug = True

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

    if is_debug:
        for k,v in qoi_manager.required_simulations.items():
            print(k)
            for k2, v2 in v.items():
                print('\t',k2,':',v2)

    sim_manager = lammps.SimulationManager()
    sim_manager.add_required_simulations(
            required_simulations = qoi_manager.get_required_simulations())
    sim_manager.structure_info = structure_info
    sim_manager.potential_info = potential_info
    if is_debug:
        print(80*"-")
        print('lammps.SimulationManager.lmps_sim_obj')
        print(80*"-")
        for k,v in sim_manager.obj_lammps_tasks.items():
            print(k,":",v)

    sim_manager.evaluate_parameters(
            param_dict = get_lewis_catlow_parameters())

    for k,v in self.obj_lammps_tasks.items():
        task_name = k
        task = v['obj']
        print(task_name)
