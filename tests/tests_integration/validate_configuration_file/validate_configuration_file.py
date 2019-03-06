import os
import importlib
import pypospack.potential
from pypospack.pyposmat.data import PyposmatConfigurationFile

if __name__ == "__main__":
    import pypospack.utils
    config_fn = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'data/Ni__eam__born_exp_fs__3.5NN/pyposmat.config.in')

    o = PyposmatConfigurationFile()
    o.read(filename=config_fn)
    o.validate()
