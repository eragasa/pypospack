import pytest

import os
import pypospack.utils

pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
configuration_filename = os.path.join(pypospack_root_dir,'data/MgO_pareto_data/pyposmat.config.in')

def test__can_import_class():
    from pypospack.pyposmat.data import PyposmatDataFile

def test____init__():
    from pypospack.pyposmat.data import PyposmatDataFile
    datafile = PyposmatDataFile()

    assert type(datafile) is PyposmatDataFile
    #assert type(datafile.names) is None
    #assert type(datafile.parameter_names) is None
    #assert type(datafile.qoi_names) is None
    #assert type(datafile.error_names) is None
    #assert type(datafile.qoi_names) is None
    #assert type(datafile.scaling_factors) is None
    #assert type(datafile.qoi_references) is None
    #assert type(datafile.scaling_factors) is None

    #assert type(datafile.df) is None
    #assert type(datafile.parameter_df) is None
    #assert type(datafile.error_df) is None
    #assert type(datafile.qoi_df) is None
    #assert type(rescaled_error_df) is None

def test____init__w_filename():
    from pypospack.pyposmat.data import PyposmatDataFile
    o = PyposmatDataFile(filename=configuration_filename)

def test__read():
    from pypospack.pyposmat.data import PyposmatDataFile
    o = PyposmatDataFile()
    o.read(filename=configuration_filename)

def dev__read():
    from pypospack.pyposmat.data import PyposmatDataFile
    o = PyposmatDataFile()
    o.read(filename=configuration_filename)
if __name__ == "__main__":
    dev_read()
