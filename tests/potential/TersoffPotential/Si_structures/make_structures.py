import pypospack.crystal as crystal
import pypospack.io.lammps as io_lammps
import pypospack.io.vasp as io_vasp

structure_dict = {
        "Ni_dia_unit":{"filename":"Ni_dia_unit.vasp",
                       "obj":None}
        }

poscar_file = io_vasp.read_poscar_file('Si_dia_unit.vasp')
lammps_file = io_lammps.LammpsStructure(poscar_file)
lammps_file.write('Si_dia_unit.lammps',atom_style='atomic')
print('done!')
