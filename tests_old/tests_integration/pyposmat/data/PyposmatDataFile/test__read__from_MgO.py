import pytest
import os
import pandas as pd
from pypospack.pyposmat.data import PyposmatDataFile

MgO_datafile = "../../../../../data/MgO_pareto_data/culled_004.out"

expected_names = ['sim_id', 'chrg_Mg', 'chrg_O', 'MgMg_A', 'MgMg_rho', 'MgMg_C', 'MgO_A', 'MgO_rho', 'MgO_C', 'OO_A', 'OO_rho', 'OO_C', 'MgO_NaCl.a0', 'MgO_NaCl.c11', 'MgO_NaCl.c12', 'MgO_NaCl.c44', 'MgO_NaCl.B', 'MgO_NaCl.G', 'MgO_NaCl.fr_a', 'MgO_NaCl.fr_c', 'MgO_NaCl.sch', 'MgO_NaCl.001s', 'MgO_NaCl.a0.err', 'MgO_NaCl.c11.err', 'MgO_NaCl.c12.err', 'MgO_NaCl.c44.err', 'MgO_NaCl.B.err', 'MgO_NaCl.G.err', 'MgO_NaCl.fr_a.err', 'MgO_NaCl.fr_c.err', 'MgO_NaCl.sch.err', 'MgO_NaCl.001s.err']

parameter_names = ['chrg_Mg', 'chrg_O', 'MgMg_A', 'MgMg_rho', 'MgMg_C', 'MgO_A', 'MgO_rho', 'MgO_C', 'OO_A', 'OO_rho', 'OO_C']

qoi_names = ['MgO_NaCl.a0', 'MgO_NaCl.c11', 'MgO_NaCl.c12', 'MgO_NaCl.c44', 'MgO_NaCl.B', 'MgO_NaCl.G', 'MgO_NaCl.fr_a', 'MgO_NaCl.fr_c', 'MgO_NaCl.sch', 'MgO_NaCl.001s']

error_names = ["{}.err".format(v) for v in qoi_names]

def test__check_harness():
    assert os.path.isfile(MgO_datafile)

def test__read():
    datafile = PyposmatDataFile()
    datafile.read(filename=MgO_datafile)

    assert type(datafile.names) is list
    assert len(expected_names) == len(datafile.names)
    for i,v in enumerate(expected_names):
        assert expected_names[i] == datafile.names[i]

    assert type(datafile.parameter_names) is list
    assert len(parameter_names) == len(datafile.parameter_names)
    for i,v in enumerate(parameter_names):
        assert parameter_names[i] == datafile.parameter_names[i]

    assert type(datafile.qoi_names) is list
    assert len(qoi_names) == len(datafile.qoi_names)
    for i,v in enumerate(qoi_names):
        assert qoi_names[i] == datafile.qoi_names[i]

    assert type(datafile.error_names) is list
    assert len(error_names) == len(datafile.error_names)
    for i,v in enumerate(error_names):
        assert error_names[i] == datafile.error_names[i]

    assert type(datafile.df) is pd.DataFrame

def test__read__wo_named_arguments():
    datafile = PyposmatDataFile()
    datafile.read(MgO_datafile)

    assert type(datafile.names) is list
    assert len(expected_names) == len(datafile.names)
    for i,v in enumerate(expected_names):
        assert expected_names[i] == datafile.names[i]

    assert type(datafile.parameter_names) is list
    assert len(parameter_names) == len(datafile.parameter_names)
    for i,v in enumerate(parameter_names):
        assert parameter_names[i] == datafile.parameter_names[i]

    assert type(datafile.qoi_names) is list
    assert len(qoi_names) == len(datafile.qoi_names)
    for i,v in enumerate(qoi_names):
        assert qoi_names[i] == datafile.qoi_names[i]

    assert type(datafile.error_names) is list
    assert len(error_names) == len(datafile.error_names)
    for i,v in enumerate(error_names):
        assert error_names[i] == datafile.error_names[i]

    assert type(datafile.df) is pd.DataFrame

if __name__ == "__main__":
    datafile = PyposmatDataFile()
    datafile.read(filename=MgO_datafile)
    print("{}:{}".format('parameter_names',datafile.parameter_names))
    print("{}:{}".format('qoi_names',datafile.qoi_names))
    print("{}:{}".format('error_names',datafile.error_names))
