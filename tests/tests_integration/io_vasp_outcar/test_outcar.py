import pypospack.io.vasp as vasp

class TestOutcar(object):

    def test_init(self):
        outcar = vasp.Outcar()

    def test_read_single(self):
        outcar = vasp.Outcar()
        outcar.read('OUTCAR.single')

        assert abs(outcar.encut 
if __name__ == '__main__':
    pass
