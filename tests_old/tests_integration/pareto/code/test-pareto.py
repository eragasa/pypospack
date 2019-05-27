from __future__ import print_function
from builtins   import range
from pareto import *
import parse

names, values = parse.read_data('../data/results_10000_a.out')
#names, values = parse.read_data('../data/20161207_results_100k.out')
vals = [[v[i] for i in range(len(v)) if names[2][i].endswith('_abserr')] for idx,p,v in values]

#timing(pareto_bruteforce)(vals)
pareto_vals = timing(pareto)(vals)

print(len(vals), len(pareto_vals))

