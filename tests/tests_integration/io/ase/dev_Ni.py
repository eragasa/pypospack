from ase.lattice.cubic import FaceCenteredCubic
import ase.build.bulk

def make_fcc_bulk(name='Ni'):
    return ase.build.bulk(name,cubic=True)

def make_fcc_111_slab(size=(2,2,3),symbol='Cu',pbc=(1,1,0)):
    direction_x = [1,-1,1]
    direction_y = [1,1,-2]
    direction_z = [1,1,1]
    atoms = FaceCenteredCubic(\
            directions=[direction_x,
                        direction_y,
                        direction_z],
            size=size,
            symbol=symbol,
            pbc=pbc)
    return atoms

def make_fcc_100_slab(size=(2,2,3),symbol='Cu',pbc=(1,1,0)):
    direction_x = [1,0,0]
    direction_y = [0,1,0]
    direction_z = [0,0,1]

    atoms = FaceCenteredCubic(\
            directions=[direction_x,
                        direction_y,
                        direction_z],
            size=size,
            symbol=symbol,
            pbc=pbc)

if __name__ == "__main__":
    make_fcc_bulk()
    import pypospack.io.vasp as vasp
    fcc_poscar = vasp.Poscar(make_fcc_bulk())
    fcc_poscar.write('Ni_fcc.vasp')
    make_fcc_111_slab()
    make_fcc_100_slab()
