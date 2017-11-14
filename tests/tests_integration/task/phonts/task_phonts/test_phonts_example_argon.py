
"""
This script tests the functionality of the pypospack.io.phonts module by 
recreating the example in the PhonTS example/argon

"""
import numpy as np
import pypospack.io.phonts as phonts
import pypospack.crystal as crystal

ar_fcc = crystal.SimulationCell()
ar_fcc.a0 = 3.9620
ar_fcc.add_atom('Ar',[0.0,0.0,0.0])
ar_fcc.add_atom('Ar',[0.5,0.5,0.0])
ar_fcc.add_atom('Ar',[0.5,0.0,0.5])
ar_fcc.add_atom('Ar',[0.0,0.5,0.5])

phonts_potential_type = ['exp-6',1]
phonts_potential_params = ['Ar','Ar',0.010531316,13.,3.85,5.625,9.625]

if __name__ == '__main__':
    phonts_ar_fcc = phonts.PhontsSimulation()




