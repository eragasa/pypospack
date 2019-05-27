import pytest
import os
from pypospack.io.vasp import Poscar

def test__import__from_pypospack_io_phonts():
    from pypospack.io.phonts import simulation_cell_to_phonts_string

def test__setup__can_read_poscar_file():
    filename = os.path.join(
            'test__simulation_cell_to_phonts_string',
            'Si_dia_unit.relax.gga.vasp')
    poscar = Poscar()
    poscar.read(filename=filename)

def test__no_charges():
    filename = os.path.join(
            'test__simulation_cell_to_phonts_string',
            'Si_dia_unit.relax.gga.vasp')
    poscar = Poscar()
    poscar.read(filename=filename)

    from pypospack.io.phonts import simulation_cell_to_phonts_string
    simulation_cell_to_phonts_string(
            simulation_cell=poscar)

if __name__ == "__main__":
    # creating tests for formatting is difficult.
    # it's better to do integration tests later to determine if
    # file created is valid.
    from pypospack.io.phonts import simulation_cell_to_phonts_string
    poscar = Poscar()
    poscar.read(filename='Si_dia_unit.relax.gga.vasp')
    cell_str= simulation_cell_to_phonts_string(
            simulation_cell=poscar)
    print(cell_str)

