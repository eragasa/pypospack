import numpy as np
import copy

position_type = ['Direct','Cartesian']

def get_reciprocal_lattice(H):
    a1 = H[0,:]
    a2 = H[1,:]
    a3 = H[2,:]

    b1 = None
    b2 = None
    b3 = None

    G = np.array([[1,0,0],[0,1,0],[0,0,1]])
    G[0,:] = b1
    G[1,:] = b2
    G[2,:] = b3

    return G
class SimulationCell(object):
    def __init__(self):
        self._coordinate_type = 'Direct'
        self._lattice = Lattice()
        self._basis = []

    @property
    def coordinate_type(self):
        """(string): coordinate type, either direct or cartesian."""
        return self._coordinate_type

    @coordinate_type.setter
    def coordinate_type(self, v):
        if v not in position_type:
            self._coordinate_type = v
        else:
            msg_out = "coordinate type must be either direct or cartesian."
            raise ValueError(msg_out)

    def add_atom(self,atom,position)

class Atom(object):
    def __init__(self):
        self._symbol = None
        self._position = None
        self._magnetic_moment = None
        self._is_interstitial = None
        self._is_vacancy = None

   
class Lattice(object):
    def __init__(self,a0=None,a1=None,a2=None,a3=None,H=None):
        self._a0 = None
        self._a1 = None
        self._a2 = None
        self._a3 = None
        self._H = None


        if (a1 is not None) and (a2 is not None) and (a3 is not None):
            for i,a in enumerate([a1,a2,a3]):
                if type(a) is list:
                   assert
            if type(a1) is list:
                self._a1 = np.array(a1)
            elif
            self._a1 = None
            self._a2
            self._a3
            H_array = [ [ a1[0], a1[1], a1[2]],
                        [ a2[0], a2[1], a2[2]],
                        [ a3[0], a2[2], a3[3]] ]
            self._H = np.array(H_array)
        else
         
