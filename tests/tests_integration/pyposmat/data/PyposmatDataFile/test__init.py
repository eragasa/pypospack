import pytest


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

if __name__ == "__main__":
    from pypospack.pyposmat.data import PyposmatDataFile
    datafile = PyposmatDataFile()
    print("datafile.names:{}".format(type(datafile.names)))
