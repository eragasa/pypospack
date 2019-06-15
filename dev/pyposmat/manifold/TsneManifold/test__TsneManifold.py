import os 

import pypospack.utils
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.manifold import Manifold
from pypospack.pyposmat.manifold import TsneManifold

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

def test____init____with_objects():
    manifold = TsneManifold(pyposmat_configuration=o_config,
                            pyposmat_data=o_data)

    assert isinstance(manifold,Manifold)

def test___init____with_paths():
    manifold = TsneManifold(pyposmat_configuration=config_fn,
                            pyposmat_data=data_fn)

    assert isinstance(manifold,Manifold)

def test__learn_manifold__qois():
    manifold = TsneManifold(pyposmat_configuration=config_fn,
                            pyposmat_data=data_fn)
    manifold.learn_manifold(names='qois')
    assert 'TSNE_1' in list(manifold.data.df.columns.values)
    assert 'TSNE_2' in list(manifold.data.df.columns.values)

def test__learn_manifold__free_parameters():
    manifold = TsneManifold(pyposmat_configuration=config_fn,
                            pyposmat_data=data_fn)
    manifold.learn_manifold(names='free_parameters')
    print(list(manifold.data.df.columns.values))
    assert 'TSNE_1' in list(manifold.data.df.columns.values)
    assert 'TSNE_2' in list(manifold.data.df.columns.values)

def test__learn_manifold__all():
    manifold = TsneManifold(pyposmat_configuration=config_fn,
                            pyposmat_data=data_fn)
    manifold.learn_manifold(names='all')
    print(list(manifold.data.df.columns.values))
    assert 'TSNE_1' in list(manifold.data.df.columns.values)
    assert 'TSNE_2' in list(manifold.data.df.columns.values)

if __name__ == "__main__":
    manifold = TsneManifold(pyposmat_configuration=config_fn,
                            pyposmat_data=data_fn)
    manifold.learn_manifold(names='all')
