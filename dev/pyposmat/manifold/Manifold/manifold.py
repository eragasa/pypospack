# -*- coding: utf-8 -*-
"""Implementation of abstract manifold class

This module provides the abstract class implementation for manifold learning
this wrapper is primarily a wrapper from the scikit-learn module.
"""
from sklearn import preprocessing

from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

from pypospack.pyposmat.manifold import PyposmatManifold

class Manifold(PyposmatManifold):
    pass
if __name__ == "__main__":
    import os
    import pypospack.utils
    # pylint: disable=invalid-name

    pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()

    config_fn = os.path.join(
        pypospack_root_dir,
        'data', 'Si__sw__data', 'pareto_optimization_unconstrained',
        'pyposmat.config.in')
        data_fn = os.path.join(
        pypospack_root_dir,
        'data', 'Si__sw__data', 'pareto_optimization_unconstrained',
        'pyposmat.kde.20.out')

    o_config = PyposmatConfigurationFile()
    o_config.read(filename=config_fn)

    o_data = PyposmatDataFile()
    o_data.read(filename=data_fn)

    o = Manifold(pyposmat_configuration=o_config,pyposmat_data=o_data)
