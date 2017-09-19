#!/bin/env python
import time
import pypospack.pyposmat as pyposmat
#import pyflamestk.potentials as potentials

"""
This script initializes the pyposmat eninge, and calculates the parameters of
interest for a set of parameters.  The parameters are print to the screen.

Required Directories:
  lmp_scripts_db  
  structure_db

Required Files:
  potential.mod
  pyposmat.config
  pyposmat.potential
  pyposmat.qoi

Output Files:
  NONE
"""
# define the set of parameters as a dictionary
param_dict = {}
param_dict['chrg_Mg'] = +2.0
param_dict['chrg_O']  = -2.0
param_dict['MgMg_A']   = 0.0 
param_dict['MgMg_rho'] = 0.5
param_dict['MgMg_C']   = 0.0
param_dict['MgO_A']    = 821.6
param_dict['MgO_rho']  = 0.3242
param_dict['MgO_C']    = 0.0
param_dict['OO_A']     = 2274.00 
param_dict['OO_rho']   = 0.1490
param_dict['OO_C']     = 27.88

start_time = time.time()

# constructor
mc_sampler = pyposmat.PyPosmatEngine()

# evaluate parameter set
names, types, values = mc_sampler.evaluate_parameter_set(param_dict)

# output results
print(names)
print(types)
print(values)
print("{} second".format(time.time()-start_time))
