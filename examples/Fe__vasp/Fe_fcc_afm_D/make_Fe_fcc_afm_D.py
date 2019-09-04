import os, shutil, subprocess
from pypospack.crystal import SimulationCell
import pypospack.io.vasp as vasp
import pypospack.crystal as crystal
import pypospack.io.slurm as minimize_init_slurm

system_name = "Fe-fccD"
# define simulation cell
cell = SimulationCell()
cell.a0 = 3.447
cell.H = [
    [1,0,0],
    [0,1,0],
    [0,0,2]
]
cell.add_atom('Fe', [0.0, 0.0, 0.00], magmom=2.22)
cell.add_atom('Fe', [0.5, 0.5, 0.00], magmom=2.22)
cell.add_atom('Fe', [0.5, 0.0, 0.25], magmom=2.22)
cell.add_atom('Fe', [0.0, 0.5, 0.25], magmom=2.22)
cell.add_atom('Fe', [0.0, 0.0, 0.50], magmom=-2.22)
cell.add_atom('Fe', [0.5, 0.5, 0.50], magmom=-2.22)
cell.add_atom('Fe', [0.5, 0.0, 0.75], magmom=-2.22)
cell.add_atom('Fe', [0.0, 0.5, 0.75], magmom=-2.22)

poscar = vasp.Poscar(cell)
poscar.write('POSCAR')
