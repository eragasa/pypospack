import pytest
from collections import OrderedDict

# The testing functions provided below are coarse-grained tests to demonstrate
# that functionality works.  More robust unit-testing of arguments generally
# isn't developed except to catch bugs in development.
#
# Imports of classes and functions should be done explicitly and within the 
# scope of each testing function.  This serves two purposes: (1) provide an
# idiomatic example of how to use the class/functions. (2) provide explicitly
# the steps reached to get expected behavior.
#
# Methods of class and functions should test return values and ensure that
# arguments which are passed in are not mutated.
# Methods of classes should not only test return values and attributes mutated
# by the method, but also insure that objects, including lists and dictionaries,
# are not mutated by the method (unless the method, explicitly does so., but also 
# expected behavior of
def test__import1():
    from pypospack.potential import EamPotential

def test_1sym____init____morse_exponential_universal():
    #<--- variables unique for the test ---------------------------------------
    symbols = ['Ni']
    func_pair='morse'
    func_density='eam_dens_exp'
    func_embedding='eam_embed_universal'

    #<--- setup of the code to conduct the test -------------------------------
    from pypospack.potential import EamPotential

    #<--- code being tested ---------------------------------------------------
    eam = EamPotential(
            symbols=symbols,
            func_pair=func_pair,
            func_density=func_density,
            func_embedding=func_embedding)

    #<--- setup testing -------------------------------------------------------
    # it isn't necessary to explicitly define the expected values of the
    # pair potentials, density functions, or embedding functions encapsulated
    # with in EamEmbeddingFunction.  Those classes should have their own suites
    # of tests developed.
    from pypospack.potential import MorsePotential
    from pypospack.potential import ExponentialDensityFunction
    from pypospack.potential import UniversalEmbeddingFunction
    pair = MorsePotential(symbols=symbols)
    dens = ExponentialDensityFunction(symbols=symbols)
    embed = UniversalEmbeddingFunction(symbols=symbols)

    p_names = ["p_{}".format(p) for p in pair.parameter_names]
    d_names = ["d_{}".format(p) for p in dens.parameter_names]
    e_names = ["e_{}".format(p) for p in embed.parameter_names]
    parameter_names = p_names + d_names + e_names
    
    #<--- setup testing attributes --------------------------------------------
    # All public attributes and properties should be tested for expected
    # behavior.  This includes all attributes and properties which are 
    # initialized to None after class constructor is called
    #<------ testing eam.obj_pair is inherited from the correct base class
    from pypospack.potential import PairPotential
    assert isinstance(eam.obj_pair,PairPotential)
    
    #<------ testing eam.obj_density is inherited from the correct base class
    from pypospack.potential import EamDensityFunction
    assert isinstance(eam.obj_density,EamDensityFunction)
    
    #<------ testing eam.obj_embedding is inherited from correct base class
    from pypospack.potential import EamEmbeddingFunction
    assert isinstance(eam.obj_embedding,EamEmbeddingFunction)
    
    #<------ testing eam.symbols
    assert type(eam.symbols) is list
    assert eam.symbols == symbols

    #<------ testing eam.parameter_names
    assert type(eam.parameter_names) is list
    assert len(eam.parameter_names) == len(parameter_names)
    for pn in parameter_names:
        pn in eam.parameter_names
    
    #<------ testing eam.parameters
    assert type(eam.parameters) is OrderedDict
    assert len(eam.parameters) == len(eam.parameter_names)
    for pn in eam.parameter_names:
        assert pn in eam.parameters
    for pn,pv in eam.parameters.items():
        assert pv is None

    #<------ testing attributes should be set to None
    assert eam.pair == None
    assert eam.density == None
    assert eam.embedding == None

def test_2sym____init____morse_exponential_universal():
    #<--- variables unique for the test ---------------------------------------
    symbols = ['Ni','Al']
    func_pair='morse'
    func_density='eam_dens_exp'
    func_embedding='eam_embed_universal'

    #<--- setup of the code to conduct the test -------------------------------
    from pypospack.potential import EamPotential

    #<--- code being tested ---------------------------------------------------
    eam = EamPotential(
            symbols=symbols,
            func_pair=func_pair,
            func_density=func_density,
            func_embedding=func_embedding)
    from pypospack.potential import EamPotential

    eam = EamPotential(
            symbols=symbols,
            func_pair='morse',
            func_density='eam_dens_exp',
            func_embedding='eam_embed_universal')

    #<--- setup testing
    from pypospack.potential import MorsePotential
    from pypospack.potential import ExponentialDensityFunction
    from pypospack.potential import UniversalEmbeddingFunction
    pair = MorsePotential(symbols=symbols)
    dens = ExponentialDensityFunction(symbols=symbols)
    embed = UniversalEmbeddingFunction(symbols=symbols)

    p_names = ["p_{}".format(p) for p in pair.parameter_names]
    d_names = ["d_{}".format(p) for p in dens.parameter_names]
    e_names = ["e_{}".format(p) for p in embed.parameter_names]
    parameter_names = p_names + d_names + e_names
    
    #<--- testing attributes
    # All public attributes and properties should be tested for expected
    # behavior.  This includes all attributes and properties which are 
    # initialized to None after class constructor is called
    #<------ testing eam.obj_pair is inherited from the correct base class
    from pypospack.potential import PairPotential
    assert isinstance(eam.obj_pair,PairPotential)
    
    #<------ testing eam.obj_density is inherited from the correct base class
    from pypospack.potential import EamDensityFunction
    assert isinstance(eam.obj_density,EamDensityFunction)
    
    #<------ testing eam.obj_embedding is inherited from correct base class
    from pypospack.potential import EamEmbeddingFunction
    assert isinstance(eam.obj_embedding,EamEmbeddingFunction)
    
    #<------ testing eam.symbols
    assert type(eam.symbols) is list
    assert eam.symbols == symbols

    #<------ testing eam.parameter_names
    assert type(eam.parameter_names) is list
    assert len(eam.parameter_names) == len(parameter_names)
    for pn in parameter_names:
        pn in eam.parameter_names
    
    #<------ testing eam.parameters
    assert type(eam.parameters) is OrderedDict
    assert len(eam.parameters) == len(eam.parameter_names)
    for pn in eam.parameter_names:
        assert pn in eam.parameters
    for pn,pv in eam.parameters.items():
        assert pv is None

    #<------ testing attributes should be set to None
    assert eam.pair == None
    assert eam.density == None
    assert eam.embedding == None
