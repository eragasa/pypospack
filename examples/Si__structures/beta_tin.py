""" make beta tin crystal structure for silicon """


from pypospack.crystal import SimulationCell
import numpy as np
class BetaTinStructure(SimulationCell):
    def __init__(self,symbols=['Pb'],a=5.8315,c=3.1814):
        assert isinstance(symbols,list)
        assert len(symbols) == 1
        SimulationCell.__init__(self)
        self.H = np.array([[a,0,0],[0,a,0],[0,0,c]])

        self.add_atom(symbol=symbols[0],position=[0.,0.,0.])
        self.add_atom(symbol=symbols[0],position=[0.5,0.0,0.75])
        self.add_atom(symbol=symbols[0],position=[0.5,0.5,0.5])
        self.add_atom(symbol=symbols[0],position=[0.0,0.5,0.25])

if __name__ == "__main__":
    from pypospack.io.vasp import Poscar
    
    # For Si, in beta-pb phase
    # 
    # 10.1103/PhysRevB.83.075119
    a = 4.81
    c_over_a = 0.552
    c = a*c_over_a

    filename = 'Si_betaPb.vasp'

    poscar = Poscar(obj_cell=BetaTinStructure(symbols=['Si'],a=a,c=c))
    poscar.write(filename=filename)

