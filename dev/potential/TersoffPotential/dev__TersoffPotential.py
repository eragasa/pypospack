import os
from collections import OrderedDict
from threebody_tersoff import get_3body_parameter_names
from threebody_tersoff import TersoffPotential

def dev__get_parameternames__1symbol():
    symbols = ['Si']
    names = TersoffPotential.twobody_parameter_names \
            + TersoffPotential.threebody_parameter_names
    parameter_names = get_3body_parameter_names(symbols=symbols,names=names)
    
    expected_count = 1 * len(names)
    print('n_expected:{}'.format(expected_count))
    print('n_actual:{}'.format(len(parameter_names)))
    #assert expected_count == len(parameter_names)


def dev__get_parameternames__2symbol():
    symbols = ['Si','Ge']

    names = TersoffPotential.twobody_parameter_names \
            + TersoffPotential.threebody_parameter_names
    parameter_names = get_3body_parameter_names(symbols=symbols,names=names)

    for i,p in enumerate(parameter_names):
        print('{:2} {}'.format(i,p))
    expected_count = 8 * len(names)
    print('n_expected:{}'.format(expected_count))
    print('n_actual:{}'.format(len(parameter_names)))
    #assert expected_count == len(parameter_names)


def dev__get_parameternames__3symbol():
    symbols = ['Si','Ge','C']
    names = TersoffPotential.twobody_parameter_names \
            + TersoffPotential.threebody_parameter_names
    parameter_names = get_3body_parameter_names(symbols=symbols,names=names)
    
    expected_count = 27 * len(names)
    print('n_expected:{}'.format(expected_count))
    print('n_actual:{}'.format(len(parameter_names)))
    #assert expected_count == len(parameter_names)

def dev__read_lammps_potential_file():
    filename = os.path.join('..','threebody','BNC.tersoff')

    symbols = ['B', 'N', 'C']

    potential = TersoffPotential(symbols=symbols)
    potential.read_lammps_potential_file(filename=filename)

    for k,v in potential.parameters.items(): print(k,v)

if __name__ == "__main__":
    dev__get_parameternames__1symbol()
    dev__get_parameternames__2symbol()
    dev__get_parameternames__3symbol()
    dev__read_lammps_potential_file()
