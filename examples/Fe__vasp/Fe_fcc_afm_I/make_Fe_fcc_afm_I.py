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

poscar = vasp.Poscar(cell)
poscar.write("POSCAR")
print(poscar.get_magmom_tag())
