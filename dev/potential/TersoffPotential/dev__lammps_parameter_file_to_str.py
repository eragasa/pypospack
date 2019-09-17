import pytest
import os
from collections import OrderedDict
from threebody_tersoff import get_3body_parameter_names
from threebody_tersoff import TersoffPotential

def dev__lammps_potential_section_to_string__BNC():
    symbols = ['B', 'N', 'C']
    potential = TersoffPotential(symbols=symbols)
    str_potential_mod = potential.lammps_potential_section_to_string()
    print('str_potential_mod:{}'.format(str_potential_mod))

def test__lammps_potential_section_to_string__BNC():
    symbols = ['B', 'N', 'C']
    potential = TersoffPotential(symbols=symbols)
    str_potential_mod = potential.lammps_potential_section_to_string()
    assert isinstance(str_potential_mod,str)

def test__lammps_potential_file__BNC():
    import pypospack.utils
    pypospack_root_path = pypospack.utils.get_pypospack_root_directory()

    symbols = ['B', 'N', 'C']
    potential = TersoffPotential(symbols=symbols)
    str_potential_mod = potential.lammps_potential_section_to_string()

if __name__ == "__main__":
    dev__lammps_potential_section_to_string__BNC()

