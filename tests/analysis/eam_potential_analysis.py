import os
from collections import OrderedDict
import pypospack.potential as potential

def get_potential_to_object_map(potential_directories):
    py_files = []
    for d in potential_directories:
        py_files += [f for f in os.listdir(d) if f.endswith('.py')]
    for f in py_files:
        print(f)

def configure_eam_potential_from_dictionary(dictionary):
    assert isinstance(dictionary,dict)
    assert dictionary['potential_type'] == 'eam'

def configure_potential_from_dictionary(dictionary):
    assert isinstance(dictionary,dict)

    if dictionary['potential_type'] == 'eam':
        return configure_eam_potential_from_dictionary(dictionary)
    else:
        pass

def get_equilibrium_density(a0,lattice):
    pass

potential_directories = ['/Users/eugeneragasa/repos/pypospack/pypospack/potentials']
get_potential_to_object_map(potential_directories)

Ni_eam__morse_exponential_universal = OrderedDict()
Ni_eam__morse_exponential_universal['symbols'] = ['Ni']
Ni_eam__morse_exponential_universal['potential_type'] = 'eam'


