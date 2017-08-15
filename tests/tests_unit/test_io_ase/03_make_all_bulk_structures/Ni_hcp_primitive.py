#!/usr/bin/env python
import numpy as np
import numpy.linalg as linalg
import ase.io
import ase.build.bulk as asebulk
import ase_module

symbol = 'Ni'
structure = 'hcp'
a = 3.508
fname_out = '{}_{}_primitive.vasp'.format(symbol,structure)

# create the orthogonal unit cell
cell = ase_module.make_structure(symbol,
        structure,
        a,
        is_cubic=False,
        is_orthorhombic=False)

# create the nonorthgonal unit cell
try:
    ase_module.write_ase_structure_as_poscar(
        ase_cell = cell, 
        fname_out=fname_out)
except FileExistsError:
    print("the file exists")
