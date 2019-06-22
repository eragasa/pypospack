from collections import OrderedDict
from pint import UnitRegistry
from pypospack.potential.eamembed_eos_zopemishin import func_zopemishin_embedding_function
from lammps_conversion import convert_to_lammps_metal_units

ureg = UnitRegistry()

potentials = OrderedDict()
potentials['density'] = OrderedDict()
potentials['density']['Ni'] = OrderedDict()
#potentials['density']['Ni']['formalism'] = func_mishin2003_density_w_cutoff
potentials['density']['Ni']['param'] = OrderedDict()
potentials['density']['Ni']['param']['r0'] = -0.3138 * ureg('nm')
potentials['density']['Ni']['param']['A0'] = 1.0 
potentials['density']['Ni']['param']['B0'] = 1.1914e4 * ureg('nanometer')
potentials['density']['Ni']['param']['C0'] = 2.0329e2 * ureg('1/nanometer**3')
potentials['density']['Ni']['param']['y'] = 1.9521e1 
potentials['density']['Ni']['param']['gamma'] = 1.6802e3 * ureg('1/nanometer')
potentials['density']['Ni']['param']['rc'] = 0.5168 * ureg('nanometer')
potentials['density']['Ni']['param']['hc'] = 0.33233 * ureg('nanometer')
potentials['density']['Ni']['param']['h0'] = 1.5 * ureg('angstrom')

potentials['pair'] = OrderedDict()
potentials['pair']['NiNi'] = OrderedDict()
#potentials['pair']['NiNi']['formalism'] = func_generalized_lj_w_cutoff 
potentials['pair']['NiNi']['param'] = OrderedDict()
potentials['pair']['NiNi']['param']['b1'] = 4.7067e-3     # no units
potentials['pair']['NiNi']['param']['b2'] = 0.15106       # no units
potentials['pair']['NiNi']['param']['r1'] = 3.8673e-3 * ureg('nanometer')
potentials['pair']['NiNi']['param']['delta'] = 3.6046e3 * ureg('eV')
potentials['pair']['NiNi']['param']['V0'] = -3.5126e3 * ureg('eV')
potentials['pair']['NiNi']['param']['rc'] = 0.5168 * ureg('nanometer')
potentials['pair']['NiNi']['param']['h'] = 0.33233 * ureg('nanometer')

potentials['embedding'] = OrderedDict()
potentials['embedding']['Ni'] = OrderedDict()
#potentials['embedding']['Ni']['formalism'] = func_zopemishin_eos
potentials['embedding']['Ni']['param'] = OrderedDict()
potentials['embedding']['Ni']['param']['a0'] = 0.3138 * 10
potentials['embedding']['Ni']['param']['B0'] = 1.1914e4 * 10
potentials['embedding']['Ni']['param']['E0'] = -4.45
potentials['embedding']['Ni']['param']['beta'] = 0.4890e-2

for k,v in potentials['density']['Ni']['param'].items():
    try:
        potentials['density']['Ni']['param'][k] = convert_to_lammps_metal_units(v)
    except AttributeError as e:
        pass
    print(k,v)
    
for k,v in potentials['pair']['NiNi']['param'].items():
    try:
        potentials['pair']['NiNi']['param'][k] = convert_to_lammps_metal_units(v)
    except AttributeError as e:
        pass
    print(k,v)
