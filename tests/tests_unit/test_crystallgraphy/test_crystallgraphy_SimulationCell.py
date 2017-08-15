import unittest
import pypospack.crystallgraphy as xtal

class TestSimulationCell_init(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_noncopy_constructor(self):
        sc = xtal.SimulationCell()

xs = xtal.SimulationCell()

if __name__ == 'main':
    unittest.main()
