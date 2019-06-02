import pytest
import os
import pypospack.crystal as crystal
import pypospack.io.vasp as vasp
import ase.build
import numpy as np

def test_Si_diamond():
    sim_cell = vasp.Poscar()
    sim_cell.read('Si_dia_unit.relax.gga.vasp')

class TestSimulationCell(object):

    def test_init_default(self):
        sim_cell = vasp.Poscar()
        assert isinstance(sim_cell,crystal.SimulationCell)

    def test_init_ase(self):
        atoms = ase.build.bulk('Cu','fcc', cubic=True)
        sim_cell = vasp.Poscar(atoms)
        assert isinstance(sim_cell,crystal.SimulationCell)

    def test_init_ase_cell(self):
        atoms = ase.build.bulk('Cu','fcc', cubic=True)
        sim_cell = vasp.Poscar(atoms)

        #check to see if we produce the correct object
        assert isinstance(sim_cell,crystal.SimulationCell)

        #check to see that lattice parameter is correct
        assert sim_cell.a0 == 1
        #check to see that the H matrix is an nd.array
        assert isinstance(sim_cell.H,np.ndarray)

        #check to see that the shape of the H matrix is correct
        assert sim_cell.H.shape == (3,3)

        #check to see that the H matrix is correct
        for i in range(3):
            for j in range(3):
                assert abs(atoms.cell[i,j]-sim_cell.H[i,j]) < 1e-6

        #check to see that the number of atoms are correct
        assert len(atoms) == len(sim_cell.atomic_basis)

        #check to see that the positions are correct
        for i,a in enumerate(atoms):
            assert sim_cell.atomic_basis[i].symbol == a.symbol
            direct_x_ase = crystal.cartesian2direct(a.position,atoms.cell)
            direct_x_pypospack = sim_cell.atomic_basis[i].position
            for j in range(3):
                assert abs(direct_x_ase[j]-direct_x_pypospack[j]) < 1e-6

    def test_init_pypospack(self):
        atoms = ase.build.bulk('Cu','fcc', cubic=True)
        sim_cell_init = vasp.Poscar(atoms)
        sim_cell_copy = vasp.Poscar(sim_cell_init)

        #check to see if we produce the correct object
        assert isinstance(sim_cell_copy,crystal.SimulationCell)

        #check to see if the lattice parameter is correct
        assert sim_cell_copy.a0 == sim_cell_init.a0

        #check to see that the H matrix is an nd.array
        assert isinstance(sim_cell_copy.H,np.ndarray)

        #check to see that the shape of the H matrix is correct
        assert sim_cell_copy.H.shape == (3,3)

        #check to see that the H matrix is correct
        for i in range(3):
            for j in range(3):
                assert abs(sim_cell_copy.H[i,j]-sim_cell_init.H[i,j]) < 1e-6

        #check to see that the number of atoms are correct
        assert len(sim_cell_init.atomic_basis) == len(sim_cell_copy.atomic_basis)

        for i,a in enumerate(sim_cell_init.atomic_basis):
            #check to see that the symbols are correct
            assert sim_cell_copy.atomic_basis[i].symbol == sim_cell_init.atomic_basis[i].symbol
            #check to see that the positions are correct
            direct_x_init = sim_cell_init.atomic_basis[i].position
            direct_x_final = sim_cell_copy.atomic_basis[i].position
            for j in range(3):
                assert abs(direct_x_init[j]-direct_x_final[j]) < 1e-6

    def test_write(self):
        atoms = ase.build.bulk('Cu','fcc', cubic=True)
        sim_cell = vasp.Poscar(atoms)
        sim_cell.write('POSCAR')
        assert os.path.exists('POSCAR')
        os.remove('POSCAR')

    def test_read(self):
        atoms = ase.build.bulk('Cu','fcc', cubic=True)
        sim_cell_init = vasp.Poscar(atoms)
        sim_cell_init.write('POSCAR')
        sim_cell_copy = vasp.Poscar()
        sim_cell_copy.read('POSCAR')

        #check to see if we produce the correct object
        assert isinstance(sim_cell_copy,crystal.SimulationCell)

        #check to see if the lattice parameter is correct
        assert sim_cell_copy.a0 == sim_cell_init.a0

        #check to see that the H matrix is an nd.array
        assert isinstance(sim_cell_copy.H,np.ndarray)

        #check to see that the shape of the H matrix is correct
        assert sim_cell_copy.H.shape == (3,3)

        #check to see that the H matrix is correct
        for i in range(3):
            for j in range(3):
                assert abs(sim_cell_copy.H[i,j]-sim_cell_init.H[i,j]) < 1e-6

        #check to see that the number of atoms are correct
        assert len(sim_cell_init.atomic_basis) == len(sim_cell_copy.atomic_basis)

        for i,a in enumerate(sim_cell_init.atomic_basis):
            #check to see that the symbols are correct
            assert sim_cell_copy.atomic_basis[i].symbol == sim_cell_init.atomic_basis[i].symbol
            #check to see that the positions are correct
            direct_x_init = sim_cell_init.atomic_basis[i].position
            direct_x_final = sim_cell_copy.atomic_basis[i].position
            for j in range(3):
                assert abs(direct_x_init[j]-direct_x_final[j]) < 1e-6


if __name__ == "__main__":
   pass 
