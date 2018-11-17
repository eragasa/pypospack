from collections import OrderedDict
import numpy as np
import pypospack.potential as potential
    
symbols = ['Mg','O']

buck = potential.BuckinghamPotential(symbols=symbols)
print('parameter_names:',buck.parameter_names)
