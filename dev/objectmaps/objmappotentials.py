import sys
import inspect

import pypospack.potential
from collections import OrderedDict

module_name = 'pypospack.potential'
cls_members = inspect.getmembers(
        pypospack.potential,
        inspect.isclass)

print(80*'-')
print('potential_subclasses')
print(80*'-')
potential_subclasses = OrderedDict()
for k,v in cls_members:
    if issubclass(v,pypospack.potential.Potential):
        potential_subclasses[k] = v
        print(v.potential_type)
exit()


print(list(potential_subclasses.keys()))
print(80*'-')
print('pair_potential_formalisms')
print(80*'-')
pair_potential_formalisms = OrderedDict()
for k,v in cls_members:
    if issubclass(v,pypospack.potential.PairPotential):
        pair_potential_formalisms[k] = v

for k,v in pair_potential_formalisms.items():
    print('{:15} {} {}'.format(k,v))

print(80*'-')
print('eam_embedding_formalisms')
print(80*'-')
eam_embedding_formalisms = OrderedDict()
for k,v in cls_members:
    if issubclass(v,pypospack.potential.EamEmbeddingFunction):
        eam_embedding_formalisms[k] = v
print(eam_embedding_formalisms)

print(80*'-')
print('eam_analytical_embedding_formalisms')
print(80*'-')
eam_analytical_embedding_formalisms = OrderedDict()
for k,v in cls_members:
    is_eam_embedding = issubclass(v,pypospack.potential.EamEmbeddingFunction)
    is_eam_eos = issubclass(v,pypospack.potential.EamEmbeddingEquationOfState)
    if is_eam_embedding and not is_eam_eos:
        eam_analytical_embedding_formalisms[k]=v
print(eam_analytical_embedding_formalisms)


print(80*'-')
print('eam_eos_embedding_formalisms')
print(80*'-')
eam_eos_embedding_formalisms = OrderedDict()
for k,v in cls_members:
    if issubclass(v,pypospack.potential.EamEmbeddingEquationOfState):
        eam_eos_embedding_formalisms[k] = v
print(eam_eos_embedding_formalisms)

print(80*'-')
print('eam_density_formalisms')
print(80*'-')
eam_density_formalisms = OrderedDict()
for k,v in cls_members:
    if issubclass(v,pypospack.potential.EamDensityFunction):
        eam_density_formalisms[k] =v
print(eam_density_formalisms)

print(80*'-')
print('threebody_potential_formalisms')
print(80*'-')
three_body_formalisms = OrderedDict()
for k,v in cls_members:
    if issubclass(v,pypospack.potential.ThreeBodyPotential):
        threebody_potential_formalisms[k] = v
print(threebody_potential_formalisms)
