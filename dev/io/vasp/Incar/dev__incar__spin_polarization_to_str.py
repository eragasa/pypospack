from pypospack.io.vasp import Incar


if __name__ == "__main__":
    incar = Incar()
    incar.ispin = 1
    print(incar._spin_polarization_to_string())

    incar.ispin = 2
    print(incar._spin_polarization_to_string())

    incar.lorbit = 10
    print(incar._spin_polarization_to_string())

    incar.lorbit = 11
    print(incar._spin_polarization_to_string())

    incar.lorbit = 12
    print(incar._spin_polarization_to_string())
