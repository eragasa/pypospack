"""
References:
  [1] Lewis, G. V., and C. R. A. Catlow. Journal of Physics C: Solid State Physics 18.6 (1985): 1149.
  [2] Martinez, Jackelyn A., et al. Computer Physics Communications 203 (2016): 201-211.
  [3] Henkelman, Graeme, et al. Physical Review B 72.11 (2005): 115437.
"""
from collections import OrderedDict

# example parameter_set
param_names = ['chrg_Mg', 'chrg_O', 
               'MgMg_A', 'MgMg_C', 'MgMg_rho',
               'OO_A',  'OO_C',  'OO_rho',
               'MgO_A', 'MgO_C', 'MgO_rho']
pot_names = ['LC+2.0','BG+2.0','BG+1.7','MA+1.7']
potentials = {}
# Lewis Catlow, from Henkelman
potentials['LC+2.0'] = [2.0, -2.0, 
                        0.0, 0.0, 0.5, 
                        22764.00, 27.88, 0.1490,   
                        821.6,  0.0, 0.3242]
# Ball Grimes, from Henkelman
potentials['BG+2.0'] = [2.0, -2.0,
                        0.0, 0.0, 0.5,
                        9547.96, 32.0, 0.21916,
                        1279.69,  0.0, 0.29969]
# Ball Grimes, from Henkelman
potentials['BG+1.7'] = [1.7, -1.7,
                        0.0, 0.0, 0.5,
                        4870, 77.0, 0.2679,
                        929.69,0.0,0.29909]
# Martinez, parameterization A
potentials['MA+1.7'] = [1.7, -1.7, 
                        0.0, 0.0, 0.5,
                        22370,17.6997, 0.28323,
                        26007,299.981, 0.21938,]
# Martinez, parameterization B
#potentials['MB+1.7'] = [1.7, -1.7,
#                          0.0, 0.0, 0.5,
#                          33174,40.6299, 0.28429,
#                          70019,288.934, 0.29122]
