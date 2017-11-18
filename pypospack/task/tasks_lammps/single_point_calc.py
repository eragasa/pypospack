from pypospack.task.lammps import LammpsSimulation

class LammpsSinglePointCalculation(LammpsSimulation):

    def _lammps_input_run_minimization(self):
        str_out = (
            '# ---- define settings\n'
            'compute eng all pe/atom\n'
            'compute eatoms all reduce sum c_eng\n'
            '# ---- run minimization\n'            
            'reset_timestep 0\n'
            'thermo 10\n'
            'thermo_style custom step pe lx ly lz xy xz yz press pxx pyy pzz pxy pxz pyz c_eatoms\n'
            'run 0\n'
            )
        return str_out

