import os
import pypospack.utils
from pypospack.pyposmat.data import (PyposmatConfigurationFile,
                                     PyposmatDataFile)
from pypospack.pyposmat.manifold import Manifold

# pypospack root directory
pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()

config_fn = os.path.join(
        pypospack_root_dir,
        'data','Si__sw__data','pareto_optimization_unconstrained',
        'pyposmat.config.in')
data_fn = os.path.join(
        pypospack_root_dir,
        'data','Si__sw__data','pareto_optimization_unconstrained',
        'pyposmat.kde.20.out')

o_config = PyposmatConfigurationFile()
o_config.read(filename=config_fn)

o_data = PyposmatDataFile()
o_data.read(filename=data_fn)

def test____init__with_objects():
    manifold = Manifold(pyposmat_configuration=o_config,
                        pyposmat_data=o_data)

def test____init__with_paths():
    manifold = Manifold(pyposmat_configuration=config_fn,
                        pyposmat_data=data_fn)

if __name__ == "__main__":
    pass


