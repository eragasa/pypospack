import os
from collections import OrderedDict
from pypospack.eamtools import SeatonSetflReader

setfl_filename = os.path.join(
        'test_EamSetflFile',
        'Ni1_Mendelev_2010.eam.fs')

setfl_reader = SeatonSetflReader(path=setfl_filename)
setfl_reader.read()
print("elements: ", setfl_reader.elements)
print("element_pairs: ", setfl_reader.element_pairs)
print("n_rho: ", setfl_reader.n_rho)
print("d_rho: ", setfl_reader.d_rho)
print("n_r: ", setfl_reader.n_r)
print("d_r: ", setfl_reader.d_r)
print("cutoff: ", setfl_reader.cutoff)
assert len(setfl_reader.embedding_function("Ni")) == setfl_reader.n_r
assert len(setfl_reader.density_function("Ni")) == setfl_reader.n_rho
assert len(setfl_reader.pair_function("NiNi")) == setfl_reader.n_r