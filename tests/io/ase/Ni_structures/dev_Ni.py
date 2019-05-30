from collections import OrderedDict
from ase.lattice.cubic import FaceCenteredCubic
from ase.lattice.cubic import BodyCenteredCubic
from ase.lattice.cubic import SimpleCubic
from ase.lattice.cubic import Diamond
from ase.lattice.hexagonal import HexagonalClosedPacked
from pypospack.io.vasp import Poscar
import numpy as np

def get_atomic_radius(symbol):
    if symbol == 'Ni':
        return  (0.5/(2**0.5))*3.524

def determine_direction_vectors_for_surface(direction=[1,1,1]):
    direction_z = np.array(direction)
    point_1 = [direction_z[0],0,0]
    point_2 = [0,direction_z[1],0]
    point_3 = [0,0,direction_z[2]]
    
    direction_x = (np.array(point_3) - np.array(point_2))
    direction_y = np.cross(
            direction_z,
            direction_x
        )
    return [ 
            direction_x.tolist(),
            direction_y.tolist(),
            direction_z.tolist()
            ]

def make_fcc_100_cell(size=(1,1,1),symbols=['Ni'],pbc=(1,1,1)):
    direction_x = [1,0,0]
    direction_y = [0,1,0]
    direction_z = [0,0,1]

    atoms = FaceCenteredCubic(\
            directions=[direction_x,direction_y,direction_z],
            size=size,
            symbol=symbols[0],
            pbc=pbc)

    return atoms

def make_fcc_111_cell(size=(1,1,1),symbols=['Ni'],pbc=(1,1,1)):
    _surface_direction = [1,1,1]

    directions = determine_direction_vectors_for_surface(
            direction=_surface_direction
        )

    size=(1,1,1)
    atoms = FaceCenteredCubic(\
            directions=directions,
            size=size,
            symbol=symbols[0],
            pbc=pbc)

    return atoms

def make_fcc_110_cell(size=(1,1,1),symbols=['Ni'],pbc=(1,1,1)):
    direction_x = [0,0,1]
    direction_y = [1,-1,0]
    direction_z = [1,1,0]

    atoms = FaceCenteredCubic(\
            directions=[direction_x,direction_y,direction_z],
            size=size,
            symbol=symbols[0],
            pbc=pbc)
    return atoms

def make_bcc_100_cell(size=(1,1,1),symbols=['Ni'],pbc=(1,1,1)):
    direction_x = [1,0,0]
    direction_y = [0,1,0]
    direction_z = [0,0,1]

    directions = [direction_x,direction_y,direction_z]
    atoms = None
    try:
        atoms = BodyCenteredCubic(\
                directions=directions,
                size=size,
                symbol=symbols[0],
                pbc=pbc)
    except ValueError as e:
        if str(e) == 'Cannot guess the bcc lattice constant of an element with crystal structure fcc.':
            r = get_atomic_radius(symbol=symbols[0])
            a0 = 4*r/(3**0.5)
            atoms = BodyCenteredCubic(\
                    directions=directions,
                    size=size,
                    symbol=symbols[0],
                    pbc=pbc,
                    latticeconstant=a0)
        else:
            raise ValueError('cannot create the bcc structure')
    return atoms

def make_sc_100_cell(size=(1,1,1),symbols=['Ni'],pbc=(1,1,1)):
    direction_x = [1,0,0]
    direction_y = [0,1,0]
    direction_z = [0,0,1]
    directions = [direction_x,direction_y,direction_z]

    atoms = None
    try:
        atoms = SimpleCubic(\
                directions=[direction_x,direction_y,direction_z],
                size=size,
                symbol=symbols[0],
                pbc=pbc)
    except ValueError as e:
        if str(e) == 'Cannot guess the sc lattice constant of an element with crystal structure fcc.':
            r = get_atomic_radius(symbol=symbols[0])
            a0 = 2*r
            atoms = SimpleCubic(\
                    directions=directions,
                    size=size,
                    symbol=symbols[0],
                    pbc=pbc,
                    latticeconstant=a0)
        else:
            raise ValueError('cannot create the simple cubic structure')
    return atoms

def make_dia_100_cell(size=(1,1,1),symbols=['Ni'],pbc=(1,1,1)):
    direction_x = [1,0,0]
    direction_y = [0,1,0]
    direction_z = [0,0,1]

    atoms = None
    try:
        atoms = Diamond(\
                directions=[direction_x,direction_y,direction_z],
                size=size,
                symbol=symbols[0],
                pbc=pbc)
    except ValueError as e:
        if str(e) == 'Cannot guess the diamond lattice constant of an element with crystal structure fcc.':
            r = get_atomic_radius(symbol=symbols[0])
            a0 = 3.567
            atoms = Diamond(\
                    directions=[direction_x,direction_y,direction_z],
                    size=size,
                    symbol=symbols[0],
                    pbc=pbc,
                    latticeconstant=a0)
        else:
            raise ValueError("cannnot create the diamond cubic structure")
    return atoms

def make_hcp_0001_cell(size=(1,1,1),symbols=['Ni'],pbc=(1,1,1)):
    direction_x = [2,-1,-1,0]
    direction_y = [0,1,-1,0]
    direction_z = [0,0,0,1]

    directions=[direction_x,direction_y,direction_z]
    atoms = None
    try:
        atoms = HexagonalClosedPacked(
            directions=[direction_x,direction_y,direction_z],
            size=size,
            symbol=symbols[0],
            pbc=pbc)
    except ValueError as e:
        if str(e) == 'Cannot guess the hcp lattice constant of an element with crystal structure fcc.':
            r = get_atomic_radius(symbol=symbols[0])
            _lattice_constants = {}
            _lattice_constants['a'] = r
            _lattice_constants['c/a'] = 1.663
            
            atoms = HexagonalClosedPacked(\
                    directions=directions,
                    size=size,
                    symbol=symbols[0],
                    pbc=pbc,
                    latticeconstant=_lattice_constants)
        else:
            raise
    return atoms

if __name__ == "__main__":
    Ni_fcc_100_unit = OrderedDict()
    Ni_fcc_100_unit['name'] = 'Ni_fcc_100'
    Ni_fcc_100_unit['symbols'] = ['Ni']
    Ni_fcc_100_unit['filename_prefix'] = 'Ni_fcc_100_unit'
    Ni_fcc_100_unit['ase'] = make_fcc_100_cell(symbols=Ni_fcc_100_unit['symbols'])
    Ni_fcc_100_unit['poscar_filename'] = '{}.vasp'.format(Ni_fcc_100_unit['filename_prefix'])
    Ni_fcc_100_unit['poscar'] = Poscar(Ni_fcc_100_unit['ase'])
    Ni_fcc_100_unit['poscar'].normalize_h_matrix()
    Ni_fcc_100_unit['poscar'].write(filename=Ni_fcc_100_unit['poscar_filename'])

    Ni_fcc_110_unit = OrderedDict()
    Ni_fcc_110_unit['name'] = 'Ni_fcc_110'
    Ni_fcc_110_unit['symbols'] = ['Ni']
    Ni_fcc_110_unit['filename_prefix'] = 'Ni_fcc_110_unit'
    Ni_fcc_110_unit['ase'] = make_fcc_110_cell(symbols=Ni_fcc_110_unit['symbols'])
    Ni_fcc_110_unit['poscar_filename'] = '{}.vasp'.format(Ni_fcc_110_unit['filename_prefix'])
    Ni_fcc_110_unit['poscar'] = Poscar(Ni_fcc_110_unit['ase'])
    Ni_fcc_100_unit['poscar'].normalize_h_matrix()
    Ni_fcc_110_unit['poscar'].write(filename=Ni_fcc_110_unit['poscar_filename'])

    Ni_fcc_111_unit = OrderedDict()
    Ni_fcc_111_unit['name'] = 'Ni_fcc_111'
    Ni_fcc_111_unit['symbols'] = ['Ni']
    Ni_fcc_111_unit['filename_prefix'] = 'Ni_fcc_111_unit'
    Ni_fcc_111_unit['ase'] = make_fcc_111_cell(symbols=Ni_fcc_111_unit['symbols'])
    Ni_fcc_111_unit['poscar_filename'] = '{}.vasp'.format(Ni_fcc_111_unit['filename_prefix'])
    Ni_fcc_111_unit['poscar'] = Poscar(Ni_fcc_111_unit['ase'])
    Ni_fcc_111_unit['poscar'].normalize_h_matrix()
    Ni_fcc_111_unit['poscar'].write(filename=Ni_fcc_111_unit['poscar_filename'])

    Ni_bcc_100_unit = OrderedDict()
    Ni_bcc_100_unit['name'] = 'Ni_bcc_100'
    Ni_bcc_100_unit['symbols'] = ['Ni']
    Ni_bcc_100_unit['filename_prefix'] = 'Ni_bcc_100_unit'
    Ni_bcc_100_unit['ase'] = make_bcc_100_cell(symbols=Ni_bcc_100_unit['symbols'])
    Ni_bcc_100_unit['poscar_filename'] = '{}.vasp'.format(Ni_bcc_100_unit['filename_prefix'])
    Ni_bcc_100_unit['poscar'] = Poscar(Ni_bcc_100_unit['ase'])
    Ni_bcc_100_unit['poscar'].normalize_h_matrix()
    Ni_bcc_100_unit['poscar'].write(filename=Ni_bcc_100_unit['poscar_filename'])

    Ni_sc_100_unit = OrderedDict()
    Ni_sc_100_unit['name'] = 'Ni_sc_100'
    Ni_sc_100_unit['symbols'] = ['Ni']
    Ni_sc_100_unit['filename_prefix'] = 'Ni_sc_100_unit'
    Ni_sc_100_unit['ase'] = make_sc_100_cell(symbols=Ni_sc_100_unit['symbols'])
    Ni_sc_100_unit['poscar_filename'] = '{}.vasp'.format(Ni_sc_100_unit['filename_prefix'])
    Ni_sc_100_unit['poscar'] = Poscar(Ni_sc_100_unit['ase'])
    Ni_sc_100_unit['poscar'].normalize_h_matrix()
    Ni_sc_100_unit['poscar'].write(filename=Ni_sc_100_unit['poscar_filename'])
    
    Ni_hcp_0001_unit = OrderedDict()
    Ni_hcp_0001_unit['name'] = 'Ni_hcp_0001_unit'
    Ni_hcp_0001_unit['symbols'] = ['Ni']
    Ni_hcp_0001_unit['filename_prefix'] = 'Ni_hcp_0001_unit'
    Ni_hcp_0001_unit['ase'] = make_hcp_0001_cell(symbols=Ni_hcp_0001_unit['symbols'])
    Ni_hcp_0001_unit['poscar_filename'] = '{}.vasp'.format(Ni_hcp_0001_unit['filename_prefix'])
    Ni_hcp_0001_unit['poscar'] = Poscar(Ni_hcp_0001_unit['ase'])
    Ni_hcp_0001_unit['poscar'].normalize_h_matrix()
    Ni_hcp_0001_unit['poscar'].write(filename=Ni_hcp_0001_unit['poscar_filename'])

    Ni_dia_100_unit = OrderedDict()
    Ni_dia_100_unit['name'] = 'Ni_dia_unit'
    Ni_dia_100_unit['symbols'] = ['Ni']
    Ni_dia_100_unit['filename_prefix'] = 'Ni_dia_100_unit'
    Ni_dia_100_unit['ase'] = make_dia_100_cell(symbols=Ni_dia_100_unit['symbols'])
    Ni_dia_100_unit['poscar_filename'] = '{}.vasp'.format(Ni_dia_100_unit['filename_prefix'])
    Ni_dia_100_unit['poscar'] = Poscar(Ni_dia_100_unit['ase'])
    Ni_dia_100_unit['poscar'].normalize_h_matrix()
    Ni_dia_100_unit['poscar'].write(filename=Ni_dia_100_unit['poscar_filename'])
    #import pypospack.io.vasp as vasp
    #fcc_poscar = vasp.Poscar(make_fcc_bulk())
    #fcc_poscar.write('Ni_fcc.vasp')
    #make_fcc_111_slab()
    #make_fcc_100_slab()
