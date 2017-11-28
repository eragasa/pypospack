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
    def __init__(self,
            task_name,
            task_directory,
            structure_filename,
            restart=False,
            fullauto=False):

        _task_type = 'lmps_min_pos'
        LammpsSimulation.__init__(self,
                task_name=task_name,
                task_directory=task_directory,
                task_type=_task_type,
                structure_filename=structure_filename,
                restart=restart,
                fullauto=fullauto)

    def _lammps_input_run_minimization(self):
        str_out = (\
            "# ---- define settings\n"
            "compute eng all pe/atom\n"
            "compute eatoms all reduce sum c_eng\n"
            "# ---- run minimization\n"
            "reset_timestep 0\n"
            "thermo 1\n"
            "thermo_style custom step pe lx ly lz xy xz yz press pxx pyy pzz pxy pxz pyz c_eatoms\n"
#            "min_style cg\n"
            "minimize 1e-20 1e-20 1000 100000\n")
        return str_out

