from collections import OrderedDict
from pint import UnitRegistry

def convert_to_lammps_metal_units(x,ureg=None):
    if ureg is None:
        ureg = UnitRegistry()
    lammps_metal_units = OrderedDict([
            ('mass','amu'),
            ('length','angstrom'),
            ('energy','electron_volt'),
            ('time','picoseconds')])

    if x.units in lammps_metal_units.values():
        return x
    str_metal_units = str(x.dimensionality)
    assert isinstance(str_metal_units,str)
    for d,u in lammps_metal_units.items():
        if d in str_metal_units:
            assert isinstance(d,str)
            assert isinstance(u,str)
            str_metal_units = str_metal_units.replace(
                    '[{}]'.format(d),
                    u)
    x.ito(ureg(str_metal_units)) 
    return x
