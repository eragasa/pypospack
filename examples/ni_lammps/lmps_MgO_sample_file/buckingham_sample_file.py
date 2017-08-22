#!/bin/env python
import pyflamestk.pyposmat
import pyflamestk.potentials
import time

# Author: Eugene Ragasa, eragasa@ufl.edu
# Date: 9/14/2016
# Description:  This file reads in a data file which contains a list of parameters, and does lammps simulations

fname_params = 'params_pub.dat'
fname_results = 'results_pub.out' 

# you should not have to modify beyond this line
 
start_time = time.time()

file_sampler = pyflamestk.pyposmat.FileParameterSampler()
file_sampler.create_lammps_simulations()
file_sampler.run(fname_params = fname_params,
                 fname_results = fname_results)

print("{} second".format(time.time()-start_time))
