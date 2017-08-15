import ase.build
import numpy as np
import pypospack.crystal as crystal
import pypospack.io.vasp as vasp

class test_poscar_init_fcc_unit(unittest.TestCase):
    
    def setUp(self):
        self.symbol = 'Cu'
        self.sg = 'fcc'
        self.a = 3.6

    def test_poscar_init_ase(unittest.TestCase):
        """ test if we can build from ase """

        obj_ppp = crystal.SimulationCell(
                ase.build.bulk(symbol,sg,a=a,cubic=True))

        poscar = vasp.Poscar(obj_ppp)

        self.assert
    def tearDown(self):
        pass


def test_poscar_init_ppp():
    """ test if we can build from pypospack.crystal """
    symbol = 'Cu'
    sg = 'fcc'
    a = 3.6

    obj_ppp_1 = crystal.SimulationCell(
            ase.build.bulk(symbol,sg,a=a,cubic=True))
    obj_ppp_2 = vasp.Poscar(obj_ppp_1)


def test_poscar_init_ase():
    """ test if we can build from ase """
    symbol = 'Cu'
    sg = 'fcc'
    a = 3.6

    obj_ppp = crystal.SimulationCell(
            ase.build.bulk(symbol,sg,a=a,cubic=True))

    poscar = vasp.Poscar(obj_ppp)

test_poscar_init_ppp()
test_poscar_init_ase()
 
