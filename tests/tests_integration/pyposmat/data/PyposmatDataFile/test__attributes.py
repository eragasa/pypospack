import pytest
import os
from pypospack.pyposmat.data import PyposmatDataFile

datafile_out_fn = "pyposmat.data.out"

def test__attribute__filename__default():
    datafile = PyposmatDataFile()
    assert datafile.filename is None

def test__attribute__filename__assignment():
    datafile = PyposmatDataFile()
    datafile.filename = datafile_out_fn
    assert datafile.filename == datafile_out_fn

def test__attribute__names__after_reading_file():
    datafile_in_fn = "../../../../../data/MgO_pareto_data/culled_004.out"
    
    datafile = PyposmatDataFile()
    datafile.read(datafile_in_fn)

    assert type(datafile.names) is list

