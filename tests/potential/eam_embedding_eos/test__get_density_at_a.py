import pytest

from collections import OrderedDict()

from pypospack.potential.eamdens_mishin2003 import func_dens_mishin2003_w_cutoff

test_args = ['a','func_density','func_density_param','lattice_type']
test_cases = OrderedDict()
test_cases['1'] = OrderedDict()
test_cases['1']['a'] = 3.52
test_cases['1']['func_density'] = func_dens_mishin2003_w_cutoff
test_cases['1']['func_density_param'] = ''
test_cases['1']['lattice_type'] = 'fcc'

from pypospack.potentia.eam_embedding_eos import get_density_at_a
@pytest.mark.parametrize(
        ",".join(test_args),
        [ for k in for test_name,test_case in test_cases
def test__get_density_at_a(
        a,func_density,func_density_param,lattice_type):
    get_density_at_a(
            a=a,
            func_density=func_density,
            func_density_param=func_density_param,
            lattice_type=lattice_type)
