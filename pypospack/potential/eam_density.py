from pypospack.potential import Potential

class EamDensityFunction(Potential):
    """implements the base class for the EAM Density Function"""

    def __init__(self,
            symbols,
            potential_type='eamdens'):

        Potential.__init__(self,
                symbols=symbols,
                potential_type=potential_type,
                is_charge=False)

        self.density_evaluations = None
