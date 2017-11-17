import pypospack.crystal as crystal
import ase.build
import numpy as np

class TestSimulationCell(object):

    def test_init_default(self):
        sim_cell = crystal.SimulationCell()
        assert isinstance(sim_cell,crystal.SimulationCell)

    def test_init_ase(self):
        atoms = ase.build.bulk('Cu','fcc', cubic=True)
        sim_cell = crystal.SimulationCell(atoms)
        assert isinstance(sim_cell,crystal.SimulationCell)

    def test_init_ase_cell(self):
        atoms = ase.build.bulk('Cu','fcc', cubic=True)
        sim_cell = crystal.SimulationCell(atoms)

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
        sim_cell_init = crystal.SimulationCell(atoms)
        sim_cell_copy = crystal.SimulationCell(sim_cell_init)

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

    def test_add_atom(self):
        sim_cell = crystal.SimulationCell()
        sim_cell.add_atom('Cu',[0,0,0])
if __name__ == "__main__":
    sim_cell = crystal.SimulationCell()
    print(sim_cell.a0)
    print(sim_cell.H)
    print(type(sim_cell.atomic_basis))
    for a in sim_cell.atomic_basis:
        print(a.symbol,a.position)

    print('add_atom')
    sim_cell.add_atom('Cu',[0,0,0])

