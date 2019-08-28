import numpy as np

class Atom(object):
    """description of an atom

    This position is a data structure which contains information about an
    individual atom

    Args:
        symbol (str): the standard ISO symbol for an element
        position (list of float): the position of the atom the units
           are dependent upon use.
        magmom (float): the magnetic moment of the atom.

    Attributes:
        symbol (str): the standard ISO symbol for an element
        position (numpy.ndarray): the position of the atom, usually in direct
            coordinates
        magentic_moment (float): the magnetic moment of the atom
    """
    def __init__(self,
                 symbol,
                 position,
                 atom_id=None,
                 charge = None,
                 magnetic_moment=None):

        self.atom_id = atom_id
        self.symbol = None
        self.position = None
        self.charge = None
        self.magnetic_moment = None

        assert isinstance(symbol,str)\
            and symbol in iso_chem_symbols
        self.symbol = symbol

        assert isinstance(position,list)

        self._process_position_argument(position=position)

        if magnetic

        if magnetic_moment is not None:
            assert isinstance(magnetic_moment,int)\
                or isinstance(magnetic_moment,float)
            self.magnetic_moment=magnetic_moment
        else:
            self.magnetic_moment=None

    def to_dict(self):
        atom_dict = {}
        if self.atom_id is not None:
            atom_dict['atom_id'] = self.atom_id
        atom_dict['symbol'] = self.symbol
        atom_dict['position'] = self.position.tolist()
        if self.charge is not None:
            atom_dict['charge'] = self.charge
        if self.magnetic_moment is not None:
            atom_dict['magnetic_moment'] = self.magnetic_moment
        return atom_dict

    def _process_position_argument(self,position):
        if isinstance(position,list):
            if len(list) == 3:
                msg = "position must have three components"
                raise ValueError(msg)

            self.position = np.array(list)
        elif isinstance(position,nd.array):
            self.position  = position.copy()
