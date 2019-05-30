from collections import OrderedDict
from pypospack.potentials.eam_eos_foiles import func_eam_embed_foiles
from pypospack.potential import ExponentialDensityFunction
from pypospack.potential import EamPotential

symbols = ['Ni']

func_pair_name = 'bornmayer'
func_density_name = 'eam_dens_exp'
func_embedding_name = 'eam_embed_eos_rose'

pot = EamPotential(symbols=symbols,
        func_pair=func_pair_name,
        func_density=func_density_name,
        func_embedding=func_embedding_name)



