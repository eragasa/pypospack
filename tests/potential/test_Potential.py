import pytest



def test__import_module():
    import pypospack.potential as potential

def test__import_class():
    from pypospack.potential import Potential

def test____import__():
    from pypospack.potential import Potential

    symbols = ['Ni']
    
    with pytest.raises(NotImplementedError):
        testpot = Potential(symbols=symbols)

if False:
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

