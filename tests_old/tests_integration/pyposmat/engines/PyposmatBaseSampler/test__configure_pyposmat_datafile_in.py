import pytest
import os

import pypospack.utils

pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
config_fn = os.path.join(pypospack_root_dir,'data/MgO_pareto_data/pyposmat.config.in')
data_in_fn = os.path.join(pypospack_root_dir,'data/MgO_pareto_data/pyposmat.kde.1.out')

def test__configure_pyposmat_datafile_in():
    from pypospack.pyposmat.engines import PyposmatBaseSampler

    sampler = PyposmatBaseSampler()
    
    assert type(sampler.data_in) is type(None)
    sampler.configure_pyposmat_datafile_in(filename=data_in_fn)

def attribute_tests(sampler):
    from pypospack.pyposmat.engines import PyposmatBaseSampler
    from pypospack.pyposmat.data import PyposmatDataFile

    assert type(sampler) is PyposmatBaseSampler
    assert type(sampler.data_in) is PyposmatDataFile
