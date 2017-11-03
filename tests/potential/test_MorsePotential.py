import pytest
from collections import OrderedDict
import numpy as np

def text_table__compare_objects_like_lists(
        obj1,
        obj2,
        name1='obj1',n
        name2='obj2'):
    
    if type(list1) not in [list]:
        raise ValueError('list1 must be a list')
    if type(list2) not is [list]:
        raise ValueError('list2 must be a list')

    # <--- cast objects into lists
    def cast_obj_into_list(obj):
        _list = None
        if isinstance(obj,list):
            _list = list(obj)
        return _list
    list1 = cast_obj_into_list(obj1)
    list2 = cast_obj_into_list(obj2)

    # <--- get the unique  elements
    all_elements = list(set(list1+list2))
    
    # <--- test
    test_elements = OrderedDict()
    for e in all_elements:
        e1 = e2 = False
        if e in list1:
            e1 = True
        if e in list2:
            e2 = True
        test_elements[e] = [e1,e2]

    # <--- print out the table
    for k,v in test_elements.items():
        print("{:20}{:20}{:20}".format



def test__import_module():
    import pypospack.potential as potential
    symbols = ['Ni']
    testpot = potential.MorsePotential(symbols)
    assert isinstance(testpot,potential.MorsePotential)

def test__import_class():
    from pypospack.potential import MorsePotential
    symbols = ['Ni']
    testpot = MorsePotential(symbols)
    assert isinstance(testpot,MorsePotential)

@pytest.fixture
def 1sym_morse(symbols=['Ni']):
    assert type(symbols) is list
    assert len(symbols) is 1

    import pypospack.potential as potential
    morse = potential.MorsePotential(symbols==symbols)
    return morse

def test__1sym____init__attr_symbols(1sym_morse):
    print('type(morse.symbols}:{}'.format(type(morse.symbols)))
    assert type(1sym_morse.symbols) is list
    assert 1sym_morse == symbols

def test__1symbol____init____():

    # <---- testing parameters
    symbols=['Ni']
    import pypospack.potential as potential
    morse = potential.MorsePotential(symbols=symbols)

    assert type(morse.symbols) is list
    assert type(morse.potential_type) is str
    assert type(morse.symbol_pairs) is list
    assert type(morse.parameter_names) is list
    assert type(morse.parameters) is OrderedDict
    assert type(morse.is_charge) is bool

    assert morse.symbols == symbols
    assert morse.potential_type == 'morse'
    assert morse.symbol_pairs == [['Ni','Ni']]
    assert morse.parameter_names == ['NiNi_D0','NiNi_a','NiNi_r0']
    
    # check the parameter dictionary
    assert len(morse.parameters) == len(morse.parameter_names)
    for p in morse.parameter_names:
        assert p in morse.parameters

def test__1symbol__evaluate__no_cutoff():
    symbols=['Ni']
    import pypospack.potential as potential
    morse = potential.MorsePotential(symbols=symbols)
    
    # <----- parameters for this test
    r_max = 11.
    N_r = 500
    r = r_max * np.linspace(1,100,N_r)/100
    
    parameters = OrderedDict()
    parameters['NiNi.D0'] = 0.001114
    parameters['NiNi.a'] = 3.429506
    parameters['NiNi.r0'] = 2.6813
    parameters['NiAl.D0'] = 0.001114
    parameters['NiAl.a'] = 3.429506
    parameters['NiAl.r0'] = 2.6813
    parameters['AlAl.D0'] = 0.001114
    parameters['AlAl.a'] = 3.429506
    parameters['AlAl.r0'] = 2.6813
    
    # <----- function being evaluated
    morse.evaluate(r,parameters,r_cut=None)

    # <----- variables to test

    # <----- test
    assert type(parameters) == OrderedDict
    assert type(potential_evaluations) == OrderedDict
    assert potential_evaluations.shape == r.shape

def test_2symbol__evaluate__no_cutoff():
    symbols = ['Ni','Al']
    import pypospack.potential as potential
    morse = potential.MorsePotential(symbols=symbols)
if False:
    def test____import__():
        from pypospack.potential import Potential
        symbols = ['Ni']
        testpot = potential.MorsePotential(symbols) 
    def test__write_lammps_potential_file():
        from pypospack.potential import Potential
        class PotentialTest(Potential):
            def __init__(self,symbols): Potential.__init__(self,symbols)
            def __init_parameter_names(self): self.parameter_names = []
            def __init_parameter_dict(self): self.parameter_dict = {}
        
        symbols = ['Ni']
        testpot = PotentialTest(symbols=symbols)

        with pytest.raises(NotImplementedError):
            testpot.write_lammps_potential_file()

    def test__write_gulp_potential_section():
        from pypospack.potential import Potential
        class PotentialTest(Potential):
            def __init_parameter_names(self): pass
            def __init_create_parameter(self): pass

        symbols = ['Ni']
        testpot = PotentialTest(symbols=symbols)

        with pytest.raises(NotImplementedError):
            testpot.write_gulp_potential_section()

    def test__lammps_potential_section_to_string():
        from pypospack.potential import Potential
        class PotentialTest(Potential):
            def __init_parameter_names(self): pass
            def __init_create_parameter(self): pass
        symbols = ['Ni']
        testpot = PotentialTest(symbols=symbols)

        with pytest.raises(NotImplementedError):
            testpot.lammps_potential_section_to_string()

    def test__gulp_potential_section_to_string():
        from pypospack.potential import Potential
        class PotentialTest(Potential):
            def __init_parameter_names(self): pass
            def __init_create_parameter(self): pass
        symbols = ['Ni']
        testpot = PotentialTest(symbols=symbols)

        with pytest.raises(NotImplementedError):
            testpot.lammps_potential_section_to_string()

    def test__phonts_potential_section_to_string():
        from pypospack.potential import Potential
        class PotentialTest(Potential):
            def __init_parameter_names(self): pass
            def __init_create_parameter(self): pass
        symbols = ['Ni']
        testpot = PotentialTest(symbols=symbols)

        with pytest.raises(NotImplementedError):
            testpot.phonts_potential_section_to_string()

