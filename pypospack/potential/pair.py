from pypospack.potential import Potential

class PairPotential(Potential):
    def __init__(self,symbols,potential_type,is_charge):
        Potential.__init__(self,
                symbols=symbols,
                potential_type=potential_type,
                is_charge=is_charge)

        self.pair_evaluations = None
 
