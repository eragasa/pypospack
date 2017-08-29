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

    def test_change_encut(self):
        incar1 = vasp.Incar()
        incar1.encut = 123
        incar1.write('INCAR.encuttest')

        incar2 = vasp.Incar(2)
        incar2.read('INCAR.encuttest')
        assert incar2.encut == incar1.encut
if __name__ == "__main__":
    filename = 'INCAR.std'

    incar = vasp.Incar()
    incar.write(filename)

    filename = 'INCAR.ibrion2'
    incar = vasp.Incar()
    incar.ibrion = 2
    incar.write(filename)
