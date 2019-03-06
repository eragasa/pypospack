import pytest
import os

from pypospack.pyposmat.data import PyposmatBadParametersFile

def test____init__wo_filename():
    f = PyposmatBadParametersFile()
    assert f.filename is None

def test__get_header_string():
    parameter_names = ['a','b','c']
    filename='test.out'
    
    f = PyposmatBadParametersFile()
    s = f.get_header_string(parameter_names=parameter_names)

    assert type(s) is str

    header_line_1 = ['sim_id'] + parameter_names + ['reason']
    header_line_2 = ['sim_id'] + len(parameter_names)*['param'] + ['reason']
    s_test = "{}\n".format(",".join(header_line_1))
    s_test += "{}\n".format(",".join(header_line_2))

    assert s_test == s

def dev__get_header_string():
    parameter_names = ['a','b','c']
    filename='test.out'
    
    f = PyposmatBadParametersFile()
    s = f.get_header_string(parameter_names=parameter_names)

    header_line_1 = ['sim_id'] + parameter_names + ['reason']
    header_line_2 = ['sim_id'] + len(parameter_names)*['param'] + ['reason']
    s_test = "{}\n".format(",".join(header_line_1))
    s_test += "{}\n".format(",".join(header_line_2))
    print(s_test)
    print(s)
    print(s_test == s)


def test__write_header_section():
    parameter_names = ['a,b,c']
    filename='test.out'
    
    if os.path.isfile(filename):
        os.remove(filename)
    
    f = PyposmatBadParametersFile()
    f.write_header_section(parameter_names=parameter_names,
                          filename=filename)

    assert os.path.isfile(filename)

    if os.path.isfile(filename):
        os.remove(filename)


if __name__ == "__main__":
    dev__get_header_string()
