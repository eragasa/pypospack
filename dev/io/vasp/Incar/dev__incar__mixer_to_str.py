from pypospack.io.vasp import Incar


if __name__ == "__main__":
    incar = Incar()
    print(incar._mixer_to_string())

    incar.amix = 1.00
    print(incar._mixer_to_string())

    incar.amix = 1.00
    incar.bmix = 0.001
    print(incar._mixer_to_string())
