""" make FCC structure for silicon """

from pypospack.crystal import SimulationCell
import numpy as np

class FaceCenteredCubicStructure(SimulationCell):
    def __init__(self,symbols=['Ni'],a=3.52):
        assert isinstance(symbols,list)
        assert len(symbols) == 1

        SimulationCell.__init__(self)
        self.H = np.array([[a,0,0],[0,a,0],[0,0,a]])

        self.add_atom(symbol=symbols[0],position=[0,0,0])
        self.add_atom(symbol=symbols[0],position=[0.5,0,0.5])
        self.add_atom(symbol=symbols[0],position=[0.5,0.5,0])
        self.add_atom(symbol=symbols[0],position=[0.0,0.5,0.5])

if __name__ == "__main__":
    from pypospack.io.vasp import Poscar

    # For Si, in fcc phase
    # Pizzagalli
    # LDA calculations
    a = 3.89

    filename = 'Si_fcc.vasp'
    poscar = Poscar(obj_cell=FaceCenteredCubicStructure(symbols=['Si'],a=a))
    poscar.write(filename=filename)

