import numpy as np
from pypospack.crystal import SimulationCell

class Diamond(SimulationCell):

    def __init__(self,symbols=['Si'],a0=5.431,cell_type='cubic'):

        SimulationCell.__init__(self)
        
        cell_initializers = {}
        cell_initializers['cubic'] = self.initialize_cubic_cell
        cell_initializers['primitive'] = self.initialize_primitive_cell

        cell_initializers[cell_type](symbol=symbols[0],a0=a0)
    
    def initialize_cubic_cell(self,symbol='Si',a0=5.431):
        self.a0 = a0
        a1 = np.array([1,0,0])
        a2 = np.array([0,1,0])
        a3 = np.array([0,0,1])
        self.H = np.stack((a1,a2,a3))
        self.add_atom(symbol=symbol,position=[0.000000, 0.000000, 0.000000])
        self.add_atom(symbol=symbol,position=[0.250000, 0.750000, 0.750000])
        self.add_atom(symbol=symbol,position=[0.500000, 0.000000, 0.500000])
        self.add_atom(symbol=symbol,position=[0.000000, 0.500000, 0.500000])
        self.add_atom(symbol=symbol,position=[0.500000, 0.500000, 0.000000])
        self.add_atom(symbol=symbol,position=[0.750000, 0.250000, 0.750000])
        self.add_atom(symbol=symbol,position=[0.750000, 0.750000, 0.250000])

    def initialize_primitive_cell(self,symbol='Si',a0=5.431):
        self.a0 = a0
        a1 = 1/2 * np.array([0,1,1])
        a2 = 1/2 * np.array([1,0,1])
        a3 = 1/2 * np.array([1,1,0])
        self.H = np.stack((a1,a2,a3))
        self.add_atom(symbol=symbol,position=[0.7500,0.75000,0.75000])
        self.add_atom(symbol=symbol,position=[0.5000,0.50000,0.50000])

if __name__ == "__main__":
    from pypospack.io.vasp import Poscar 
    
    o = Poscar(obj_cell=Diamond(symbols=['Si'],a0=5.431,cell_type='cubic'))
    o.write(filename="Si_dia_unit.vasp")

    o = Poscar(obj_cell=Diamond(symbols=['Si'],a0=5.431,cell_type='primitive'))
    o.write(filename='Si_dia_prim.vasp')

    o = Diamond(symbols=['Si'],a0=5.431,cell_type='primitive')
    print(o.b1)
    print(o.b2)
    print(o.b3)

