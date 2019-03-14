import ase.build
import numpy as np
from pypospack.crystal import SimulationCell

class DiamondStructure(SimulationCell):

    def __init__(self,symbols=['Si'],a=5.43):

        assert isinstance(symbols,list)
        assert len(symbols) == 1

        SimulationCell.__init__(self)
        self.H = np.array([[a,0,0],[0,a,0],[0,0,a]])
        self.add_atom(symbol=symbols[0],position=[0.0,0.0,0.0])
        self.add_atom(symbol=symbols[0],position=[0.25,0.25,0.25])
        self.add_atom(symbol=symbols[0],position=[0.00,0.50,0.50])
        self.add_atom(symbol=symbols[0],position=[0.25,0.75,0.75])
        self.add_atom(symbol=symbols[0],position=[0.50,0.00,0.50])
        self.add_atom(symbol=symbols[0],position=[0.75,0.25,0.75])
        self.add_atom(symbol=symbols[0],position=[0.50,0.50,0.00])
        self.add_atom(symbol=symbols[0],position=[0.75,0.75,0.25])

if __name__ == "__main__":
    from pypospack.io.vasp import Poscar
    filename = "Si_dia_unit.vasp"

    cell = DiamondStructure(symbols=['Si'],a=5.43)
    poscar = Poscar(obj_cell=DiamondStructure(symbols=['Si'],
                                              a=5.43)
    )
    poscar.write(filename=filename)
