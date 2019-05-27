import pytest
import os
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile

_data_fn = os.path.join('rank_0','pyposmat.results.out')

def test____init__no_args():
    o_data = PyposmatDataFile()
    o_data.read(filename=_data_fn)

def test____init____with_filename():
    o_data = PyposmatDataFile(filename=_data_fn)

    assert o_data.filename == _data_fn

def test__read__with_filename():
    o_data = PyposmatDataFile()
    o_data.read(filename=_data_fn)

    # check number of columns
    assert type(o_data.n_samples) is int
if __name__ == "__main__":
    o_data = PyposmatDataFile()
    o_data.read(filename=_data_fn)
    (n_rows,n_cols) = o_data.df.shape
    print(o_data.df['sim_id'].max())
    print(o_data.n_samples)
