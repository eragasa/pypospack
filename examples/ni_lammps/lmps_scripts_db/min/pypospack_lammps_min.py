import pypospack.task.vasp as vasptasks
import pypospack.task.lammps as lammpstasks

fn_poscar = None

class Potential(object):
    def __init__(self):
        pass

    def clone(self):
        pass

potential = Potential()

o_testmin = lammpstasks.LammpsMinimizeStructure(fn_poscar,potential)
o_testmin.write_lammps_input_file()
