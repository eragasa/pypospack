from pypospack.task.lammps import LammpsSimulation

class LammpsPositionMinimization(LammpsSimulation):
    """ Class for LAMMPS position minimization

    This class sets up, runs and processes the relaxation of atomic positions
    to find the lowest energy structure in the local basin of attraction.  The
    simulation cell is frozen.

    Args:
        task_name (str): the name of the task
        task_directory (str): the directory of the task
    """
    def __init__(self,task_name,task_directory):
        LammpsSimulation.__init__(self,task_name,task_directory)

    def ready(self,task_dict):
        LammpsSimulation.ready(self,task_dict)
    
    def run(self, param_dict):
        LammpsSimulation.run(self,param_dict)

    def modify_structure(self):
        # i need to figure out how to modify this structure pre-simulation
        a1 = None
        a2 = None
        a3 = None
        for k,v in self.task_dict.items():
            structure_name = k.strip().split('.')[0]
            simulation_type = k.strip().split('.')[1]
            simulation_name = "{}.{}".format(structure_name,simulation_type)
            variable_name = k.strip().split('.')[2]
            if variable_name == 'xx':
                a1 = v
            elif variable_name == 'yy':
                a2 = v
            elif variable_name == 'zz':
                a3 = v
            else:
                pass

        cond1 = a1 == a2
        cond2 = a2 == a3
        cond3 = a1 == a3
        if cond1 and cond2 and cond3:
            self.lammps_structure.a0 == a1

    def lammps_input_file_to_string(self):
        str_out = "".join([\
                self._lammps_input_initialization_section(),
                self._lammps_input_create_atoms(),
                self._lammps_input_define_potential(),
                self._lammps_input_position_minimization(),
                self._lammps_input_out_section()])
        return(str_out)

    def _lammps_input_position_minimization(self):
        str_out = (\
            "# ---- define settings\n"
            "compute eng all pe/atom\n"
            "compute eatoms all reduce sum c_eng\n"
            "# ---- run minimization\n"
            "reset_timestep 0\n"
            "thermo 1\n"
            "thermo_style custom step pe lx ly lz xy xz yz press pxx pyy pzz pxy pxz pyz c_eatoms\n"
            "min_style cg\n"
            "minimize 1e-20 1e-20 1000 100000\n")
        return str_out
