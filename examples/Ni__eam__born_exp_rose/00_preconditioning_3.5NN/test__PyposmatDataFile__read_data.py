import pytest

import os
import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataFile

config_fn =  os.path.join(
        pypospack.utils.get_pypospack_root_directory(),
        'examples/Ni__eam__born_exp_rose/preconditioning_3.5NN/data/pyposmat.config.in'
        )

data_fn = os.path.join(
        pypospack.utils.get_pypospack_root_directory(),
        'examples/Ni__eam__born_exp_rose/preconditioning_3.5NN/rank_0/pyposmat.results.out'
        )

def test__read_data():
    o_data = PyposmatDataFile()
    o_data.read(filename=data_fn)

if __name__ == "__main__":
    o_data = PyposmatDataFile()
    o_data.read(filename=data_fn)
