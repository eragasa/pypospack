import pytest
import os
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatBadParametersFile

def get_testing_set():
    testing_set = OrderedDict()
    testing_set['parameter_names'] = ['a','b','c']
    testing_set['filename'] = 'test.out'
    testing_set['parameters'] = OrderedDict([
        ('a',1.),
        ('b',2.),
        ('c',3.)
    ])
    return testing_set

def test__get_header_string():
    f = PyposmatBadParametersFile()
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


def dev__write_header_section():
    print(80*'-')
    print('{:^80}'.format('write_header_section'))
    print(80*'-')

    testing_set = get_testing_set()
    
    if os.path.isfile(testing_set['filename']):
        os.remove(testing_set['filename'])

    f = PyposmatBadParametersFile()
    f.write_header_section(parameter_names=testing_set['parameter_names'],
                           filename=testing_set['filename'])

    with open(testing_set['filename'],'r') as f:
        print(f.read())

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

def dev__write_simulation_exception():
    import importlib
    
    testing_set = get_testing_set()

    o = PyposmatBadParametersFile()
    o.write_header_section(parameter_names=testing_set['parameter_names'],
                           filename=testing_set['filename'])
    
    module_name = 'pypospack.exceptions'
    module = importlib.import_module(module_name)
    
    class_names = ['LammpsSimulationError']
    for class_name in class_names:
        sim_id = class_name
        m = "message"
        exception = getattr(module,class_name)(m,parameters=testing_set['parameters'])
        o.write_simulation_exception(sim_id=sim_id,exception=exception)



if __name__ == "__main__":
    dev__get_header_string()
    dev__write_header_section()
    dev__write_simulation_exception()
