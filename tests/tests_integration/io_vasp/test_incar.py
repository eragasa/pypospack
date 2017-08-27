import pypospack.io.vasp as vasp

if __name__ == "__main__":
    filename = 'INCAR.std'

    incar = vasp.Incar()
    incar.write(filename)

    filename = 'INCAR.ibrion2'
    incar = vasp.Incar()
    incar.ibrion = 2
    incar.write(filename)
