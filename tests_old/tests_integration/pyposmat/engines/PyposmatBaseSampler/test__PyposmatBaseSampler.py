import pytest

import os
from pathlib import Path

import pypospack.utils
from pypospack.pyposmat.data import PyposmatConfigurationFile

pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
config_fn = os.path.join(pypospack_root_dir,'data/MgO_pareto_data/pyposmat.config.in')

def test__MgO_configuration_file_exists():
    assert os.path.isdir(pypospack_root_dir)
    assert os.path.isfile(config_fn)

if __name__ == "__main__":
    from pypospack.pyposmat.engines import PyposmatBaseSampler
    sampler = PyposmatBaseSampler()

