from pypospack.io.vasp import Poscar
from pypospack.crystal import make_super_cell

_bulk_filename =  "Ni_fcc_111_unit.gga.relaxed.vasp"
_bulk_poscar = Poscar()
_bulk_poscar.read(_bulk_filename)
_isf_filename = "Ni_fcc_isf.vasp"
print('a1',_bulk_poscar.a1)
print('a1',_bulk_poscar.a2)
print('a1',_bulk_poscar.a3)

_isf_sc=[1,1,5]
_isf_sc=Poscar(
        make_super_cell(_bulk_poscar,_isf_sc)
    )
import numpy as np
_isf_sc.H[2,2] =_isf_sc.H[2,2] * 0.933333
for a in _isf_sc.atomic_basis:
    a.position[2] = a.position[2] / 0.933333
_isf_sc.write(_isf_filename)    
