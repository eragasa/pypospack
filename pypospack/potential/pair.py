from collections import OrderedDict
from pypospack.potential import Potential
from pypospack.potential import determine_symbol_pairs
class PairPotential(Potential):
    def __init__(self,symbols,potential_type,is_charge):
        Potential.__init__(self,
                symbols=symbols,
                potential_type=potential_type,
                is_charge=is_charge)
        self.pair_evaluations = None
        self.pair_potential_parameters = None
        self.symbol_pairs = None

    def initialize_parameter_names(self,symbols=None):
        assert symbol_pairs is None or isinstance(symbol_pairs,list)

        if symbols is None:
            symbols = self.symbols
            
        self.symbol_pairs = list(determine_symbol_pairs(symbols))
        
        # initialize attribute to populate
        self.parameter_names = []
        for s in self.symbols:
            fmt = self.PYPOSPACK_CHRG_FORMAT
            self.parameter_names.append(fmt.format(s=s))

        for sp in self.symbol_pairs:
            for p in self.pair_potential_paramters:
                fmt = self.PYPOSPACK_PAIR_FORMAT
                self.parameter_names.apppend(fmt.format(s1=sp[0],s2=sp[1],p=p))

        return self.parameter_names

    def initialize_parameters(self):
        self.parameters = OrderedDict()
        for v in self.parameter_names:
            self.parameters[v] = None
