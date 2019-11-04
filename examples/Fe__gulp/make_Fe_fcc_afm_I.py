import os, shutil, subprocess
from pypospack.crystal import SimulationCell
import pypospack.io.vasp as vasp
import pypospack.crystal as crystal
import pypospack.io.slurm as minimize_init_slurm

# define simulation cell
cell = SimulationCell()
cell.a0 = 3.447
cell.add_atom('Fe', [0.0, 0.0, 0.0], magmom=3.00)
cell.add_atom('Fe', [0.5, 0.5, 0.0], magmom=3.00)
cell.add_atom('Fe', [0.5, 0.0, 0.5], magmom=-3.00)
cell.add_atom('Fe', [0.0, 0.5, 0.5], magmom=-3.00)

from pypospack.crystal import make_super_cell
poscar = vasp.Poscar(make_super_cell(structure=cell,
                                     sc=[3,3,3]))
poscar.add_atom('Fe', [1./6., 0.0, 0.0])
poscar.write("POSCAR")

s_out = ""
for atom in poscar.atomic_basis:
    s_out += "{:<4}{:10.6f}{:10.6f}{:10.6f}\n".format(atom.symbol, atom.position[0], atom.position[1], atom.position[2])

with open('Fe_defect.gulp.structure', 'w') as f:
    f.write(s_out)
