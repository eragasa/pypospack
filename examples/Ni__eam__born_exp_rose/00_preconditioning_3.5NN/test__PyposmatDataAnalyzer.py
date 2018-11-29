import pytest, os
import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataAnalyzer
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile

pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
pyposmat_config_fn = os.path.join(
        pypospack_root_dir,
        'examples/Ni__eam__born_exp_rose/preconditioning_3.5NN/data/pyposmat.config.in')
pyposmat_data_fn = os.path.join(
        pypospack_root_dir,
        'examples/Ni__eam__born_exp_rose/preconditioning_3.5NN/data/pyposmat.results.2.out')
pyposmat_kde_fn = os.path.join(
        'pyposmat.kde.out')

def test__init():
    o = PyposmatDataAnalyzer()

def test__init__with_filenames():
    o = PyposmatDataAnalyzer(
            fn_config=pyposmat_config_fn,
            fn_data=pyposmat_data_fn)

    assert type(o.configuration) is PyposmatConfigurationFile
    assert type(o.data) is PyposmatDataFile

def test__read_configuration():
    o = PyposmatDataAnalyzer()
    o.read_configuration(filename=pyposmat_config_fn)

    assert type(o.configuration) is PyposmatConfigurationFile

def test__read_data():
    o = PyposmatDataAnalyzer()
    o.read_data(filename=pyposmat_data_fn)

    assert type(o.data) is PyposmatDataFile

def test__write_kde_file():
    o = PyposmatDataAnalyzer()
    o.read_configuration(filename=pyposmat_config_fn)
    o.read_data(filename=pyposmat_data_fn)
    o.write_kde_file(filename=pyposmat_kde_fn)

    assert os.path.isfile(pyposmat_kde_fn)

    # cleanup
    os.remove(pyposmat_kde_fn)

if __name__ == "__main__":
    o = PyposmatDataAnalyzer()
    o.read_configuration(filename=pyposmat_config_fn)
    o.read_data(filename=pyposmat_data_fn)
    o.write_kde_file(filename=pyposmat_kde_fn)
