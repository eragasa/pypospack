__author__ = "Eugene J. Ragasa, Xue Xiong"

from ase.lattice.hexagonal import HexagonalFactory

class HexagonalSiO2Factory(HexagonalFactory):
    """
    Args:
        directions
        size
        symbol
        lattice_constant
    element_basis(tuple of int): 1 is Si, 0 is O
    """
    bravais_basis = [ [0.476305, 0.476305, 0.000000], #Si
                      [0.523695, 0.000000, 0.333333], #Si
                      [0.000000, 0.523695, 0.666667], #Si
                      [0.585146, 0.743906, 0.205457], #O
                      [0.158760, 0.414854, 0.872123], #O
                      [0.256094, 0.841240, 0.538790], #O
                      [0.743906, 0.585146, 0.794543], #O
                      [0.841240, 0.256094, 0.461210], #O
                      [0.414854, 0.158760, 0.127877]] #O
    element_basis = (0,0,0,1,1,1,1,1,1)

SiO2_a_quartz = HexagonalSiO2Factory()

if __name__ == "__main__":
    # directions for hexagonal crystals must be specificed by the 
    # miller indices
    # a1, a2, a3, and z
    # a1 = [a11,a12,a13]
    # a2 = [a21,a22,a23]
    # a3 = [a31,a32,a33]
    # where sum(a11,a12,a13) == 0
    #       sum(a21,a22,a23) == 0
    #       sum(a31,a32,a33) == 0
    _directions = ( [2,-1,-1,0], # [1,0,0] -> x-axis equivalent in miller vector
                    [0,1,-1,0],  # [0,1,0] -> y-axis equivalent in miller vector
                    [0,0,0,1] )  # [0,0,1] -> z-axis equivalent in miller vector
    _latticeconstant = {} 
    _latticeconstant['a'] = 5.022 # Angs
    _latticeconstant['b'] = 5.022 # Angs
    _latticeconstant['c'] = 5.511 # Angs
    _latticeconstant['alpha'] = 90 #deg
    _latticeconstant['beta'] = 90 # degrees
    _latticeconstant['gamma'] = 120 # degrees

    # make ASE structure
    from ase import Atoms, Atom
    SiO2_aq = SiO2_a_quartz(directions=_directions,
                            size=(1,1,1),
                            symbol=('Si','O'),
                            latticeconstant=_latticeconstant)

    # convert to pypospack.crystal.SimulationCell
    from pypospack.crystal import SimulationCell
    SiO2_pyp = SimulationCell(SiO2_aq)

    # write the poscar file
    from pypospack.io.vasp import Poscar
    SiO2_poscar = Poscar(SiO2_pyp)
    SiO2_poscar.write('SiO2_aq_cubic.vasp')

    # write lammps structure
    from pypospack.io.lammps import LammpsStructure
    SiO2_lammps = LammpsStructure(SiO2_pyp)
    SiO2_lammps.write(filename='SiO2_aq_cubic.lammps',
            symbol_list=['Si','O'],
            atom_style='charge') #or 'atomic'
