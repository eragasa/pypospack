import pytest
import os,copy
from collections import OrderedDict

potential_configuration_MgO = OrderedDict()
potential_configuration_MgO['potential_type'] = 'buckingham'
potential_configuration_MgO['symbols'] = ['Mg','O']
potential_configuration_MgO['parameter_names'] = [
        'chrg_Mg', 'chrg_O',
        'MgMg_A', 'MgMg_rho', 'MgMg_C',
        'MgO_A','MgO_rho', 'MgO_C',
        'OO_A','OO_rho', 'OO_C']

parameter_definition = OrderedDict()
parameter_definition['type'] = 'distribution'
parameter_definition['chrg_Mg'] = OrderedDict()
parameter_definition['chrg_Mg']['distribution'] = 'uniform'
parameter_definition['chrg_Mg']['low'] = 1.5
parameter_definition['chrg_Mg']['high'] = 2.5
parameter_definition['chrg_O'] = OrderedDict()
parameter_definition['chrg_O'] = {}
parameter_definition['chrg_O']['equals'] = '-chrg_Mg'
parameter_definition['MgMg_A'] = OrderedDict()
parameter_definition['MgMg_A']['equals'] = 0.0
parameter_definition['MgMg_rho'] = OrderedDict()
parameter_definition['MgMg_rho']['equals'] = 0.5
parameter_definition['MgMg_C'] = OrderedDict()
parameter_definition['MgMg_C']['equals'] = 0.0
parameter_definition['MgO_A'] = OrderedDict()
parameter_definition['MgO_A']['distribution'] = 'uniform'
parameter_definition['MgO_A']['low'] = 800.0
parameter_definition['MgO_A']['high'] = 1300.0
parameter_definition['MgO_rho'] = OrderedDict()
parameter_definition['MgO_rho']['distribution'] = 'uniform'
parameter_definition['MgO_rho']['low'] = 0.2900
parameter_definition['MgO_rho']['high'] = 0.3300
parameter_definition['MgO_C'] = OrderedDict()
parameter_definition['MgO_C']['equals'] = 0.00
parameter_definition['OO_A'] = OrderedDict()
parameter_definition['OO_A']['distribution'] = 'uniform'
parameter_definition['OO_A']['low'] = 500.00
parameter_definition['OO_A']['high'] = 25000.00
parameter_definition['OO_rho'] = OrderedDict()
parameter_definition['OO_rho']['distribution'] = 'uniform'
parameter_definition['OO_rho']['low'] = 0.1000
parameter_definition['OO_rho']['high'] = 0.4000
parameter_definition['OO_C'] = OrderedDict()
parameter_definition['OO_C']['distribution'] = 'uniform'
parameter_definition['OO_C']['low'] = 25.0000
parameter_definition['OO_C']['high'] = 77.0000

potential_configuration_MgO = OrderedDict()
potential_configuration_MgO['potential_type'] = 'buckingham'
potential_configuration_MgO['symbols'] = ['Mg','O']
potential_configuration_MgO['parameter_names'] = [
        'chrg_Mg', 'chrg_O',
        'MgMg_A', 'MgMg_rho', 'MgMg_C',
        'MgO_A','MgO_rho', 'MgO_C',
        'OO_A','OO_rho', 'OO_C']
potential_configuration_MgO['parameter_definitions'] \
        = copy.deepcopy(parameter_definition) 
potential_configuration_MgO['free_parameter_names'] \
        = ['chrg_Mg','MgO_A','MgO_rho','OO_A','OO_rho','OO_C']

def test__import__from_pypospack_potential():
    from pypospack.potential import PotentialInformation
     
def test____init__():
    from pypospack.potential import PotentialInformation
    potinfo = PotentialInformation()

def test__set_attribute__potential_type():
    from pypospack.potential import PotentialInformation
    potinfo = PotentialInformation()

    potential_type = potential_configuration_MgO['potential_type']

    potinfo.potential_type = potential_type

    assert potinfo.potential_type == potential_type
    assert potinfo.eam_pair is None
    assert potinfo.eam_embedding is None
    assert potinfo.eam_density is None
    
def test__set_attribute__potential_type():
    from pypospack.potential import PotentialInformation
    potinfo = PotentialInformation()
   
    symbols = potential_configuration_MgO['symbols']

    potinfo.symbols = symbols
    
    assert potinfo.symbols == symbols

def test__set_attribute__parameter_definitions():
    from pypospack.potential import PotentialInformation
    potinfo = PotentialInformation()
    
    
    parameter_names = potential_configuration_MgO['parameter_names']
    parameter_definitions = potential_configuration_MgO['parameter_definitions']
    
    potinfo.parameter_names = parameter_names
    potinfo.parameter_definitions = parameter_definitions

    free_parameter_names = potential_configuration_MgO['free_parameter_names']
    assert potinfo.free_parameter_names == free_parameter_names

def test__write():
    yaml_filename_out = os.path.join(
            'pypospack.potential.yaml')

    from pypospack.potential import PotentialInformation
    potinfo = PotentialInformation()
    potinfo.write(
            filename=yaml_filename_out,
            potential_configuration=potential_configuration_MgO)

def test__read():
    yaml_filename_in = os.path.join(
            'test_PotentialConfiguration',
            'pypospack.potential.buckingham.yaml')

    from pypospack.potential import PotentialInformation
    potinfo = PotentialInformation()
    potinfo.read(filename=yaml_filename_in)

    symbols = potential_configuration_MgO['symbols']
    potential_type = potential_configuration_MgO['potential_type']
    parameter_names = potential_configuration_MgO['parameter_names'] 

    assert potinfo.symbols == symbols
    assert potinfo.potential_type == potential_type
    assert potinfo.parameter_names == parameter_names
if __name__ == '__main__':
    from pypospack.potential import PotentialInformation
    
    # create buckingham yaml file for testing
    buckingham_yaml_filename = os.path.join(
            'test_PotentialConfiguration',
            'pypospack.potential.buckingham.yaml')
    potinfo = PotentialInformation()
    potinfo.write(
            filename=buckingham_yaml_filename,
            potential_configuration=potential_configuration_MgO)

    # create eam yaml file for testing

    # check to see if can read what we wrote
    #copy_pot_info = potfit.PotentialInformation()
    #copy_pot_info.read(fname='pypospack.potential.yaml')

    #pot_info = potfit.PotentialInformation()
    #pot_info.read(fname='pypospack.buckingham.yaml')

