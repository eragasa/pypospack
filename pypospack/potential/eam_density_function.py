from pypospack.potential import Potential

class EamDensityFunction(Potential):

    potential_type = 'eam_density_base'
    def __init__(self,
            symbols,
            potential_type='eamdens'):

        Potential.__init__(self,
                symbols=symbols,
                potential_type=potential_type,
                is_charge=False)

        self.density_evaluations = None

# from pypospack.potentials.eam_dens_exp import ExponentialDensityFunction
