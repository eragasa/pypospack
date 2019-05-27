#!/usr/bin/env python
import ase.io
import ase.build.bulk as asebulk

def make_structure(symbol,structure,a,cubic=True):
    """ creates an atomic structure
    Args:
        symbol(str): the ISO chemical symbol for atom type
        structure(str): structure type
        a(float): lattice parameter for the material
        cubic(bool): true
    Returns:
        ase.atom.Atoms
    """
    sim_cell = asebulk(symbol, structure, a, cubic=cubic)
    return sim_cell

def write_ase_structure_as_poscar(ase_cell, fname_out='POSCAR'):
    """ write the ase structure
    Args:
       ase_cell(ase.atom.Atoms): the simulation cell
       fname(str): file name to be written to (default: 'POSCAR')
    """
    ase.io.write(
           filename=fname_out,
           images=ase_cell,
           format='vasp')

if __name__ == "__main__":
    symbol = 'Ni'
    structure = 'fcc'
    a = 3.508

    # create the orthogonal unit cell
    cell = make_structure(symbol,
            structure,
            a,
            cubic=True)
    # create the nonorthgonal unit cell
    write_ase_structure_as_poscar(
            ase_cell = cell, 
            fname_out='Ni_fcc_unit.vasp')

    
