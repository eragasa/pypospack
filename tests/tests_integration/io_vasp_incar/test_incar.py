import pypospack.io.vasp as vasp

class TestIncar(object):

    def test_init(self):
        incar = vasp.Incar()

    def test_write(self):
        incar = vasp.Incar()
        incar.write('INCAR.std')

    def test_read(self):
        incar = vasp.Incar()
        incar.read('INCAR.std')

    def test_write_w_ibrion_2(self):
        incar = vasp.Incar()
        incar.ibrion = 2
        incar.isif = 3
        incar.write('INCAR.ibrion2')

if __name__ == "__main__":
    filename = 'INCAR.std'

    incar = vasp.Incar()
    incar.write(filename)

    filename = 'INCAR.ibrion2'
    incar = vasp.Incar()
    incar.ibrion = 2
    incar.write(filename)
