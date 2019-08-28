from pint import UnitRegistry
from lammps_conversion import convert_to_lammps_metal_units

if __name__ == "__main__":
    ureg = UnitRegistry()
    a =   -0.3138 * ureg('nm')
    
    convert_to_lammps_metal_units(a)
    print(a)

