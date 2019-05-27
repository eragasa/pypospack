from pypospack.io.vasp import Poscar
from pypospack.crystal import make_super_cell

_bulk_filename =  "Ni_fcc_111_unit.gga.relaxed.vasp"
_bulk_poscar = Poscar()
_bulk_poscar.read(_bulk_filename)
_surface_filename = "Ni_fcc_111_surface.vasp"
print(_bulk_poscar.a1)
print(_bulk_poscar.a2)
print(_bulk_poscar.a3)

_surface_filename = "Ni_fcc_111_surf.vasp"
_surface_sc=[1,1,6]
_surface_sc=Poscar(
        make_super_cell(_bulk_poscar,_surface_sc)
    )
_surface_slab_pct_of_z = 0.5
_surface_slab_thickness = _surface_sc.a3*_surface_slab_pct_of_z
print("surface_slab_thickness={}".format(_surface_slab_thickness))

atoms_to_remove = []
for a in _surface_sc.atomic_basis:
    print(a.position,a.position[2] < _surface_slab_pct_of_z)
    if not a.position[2] < _surface_slab_pct_of_z:
        atoms_to_remove.append([a.symbol,list(a.position)])

for a in atoms_to_remove:
    symbol = a[0]
    position = a[1]
    print('removing {} atom @ {}'.format(
            str(a[0]),
            str(a[1]))
        )
    _surface_sc.remove_atom(symbol=symbol,position=position)

for a in _surface_sc.atomic_basis:
   print(a.symbol,a.position)

_surface_sc.write(_surface_filename)
