import os
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatDataFile

datafile_out_fn = "pyposmat.data.out"


def cleanup_test():
    if os.path.isfile(datafile_out_fn):
        os.remove(datafile_out_fn)

def test__get_header_string():
    
    parameter_names = ['param{}'.format(i+1) for i in range(3)]
    qoi_names = ['qoi{}'.format(i+1) for i in range(5)]
    error_names = ['err{}'.format(i+1) for i in range(5)]
    
    datafile = PyposmatDataFile()
    s = datafile.get_header_string(
            parameter_names = parameter_names,
            qoi_names = qoi_names,
            error_names = error_names)

    assert type(s) is str
    
    # check assignment
    assert type(datafile.parameter_names) is list
    assert len(datafile.parameter_names) == len(parameter_names)
    for i,v in enumerate(parameter_names):
        assert datafile.parameter_names[i] == v

    assert type(datafile.qoi_names) is list
    assert len(datafile.qoi_names) == len(qoi_names)
    for i,v in enumerate(qoi_names):
        assert datafile.qoi_names[i] == v

    assert type(datafile.error_names) is list
    assert len(datafile.error_names) == len(error_names)
    for i,v in enumerate(error_names):
        assert datafile.error_names[i] == v

    # check string
    lines = s.split("\n")
    line_1 = lines[0].strip().split(",")
    line_2 = lines[1].strip().split(",")

    # check line 1
    assert 'sim_id' in line_1
    assert 'cluster_id' not in line_1
    for v in parameter_names: assert v in line_1
    for v in qoi_names: assert v in line_1
    for v in error_names: assert v in line_1

    # check line 2
    assert line_2.count('sim_id') == 1
    assert line_2.count('cluster_id') == 0
    assert line_2.count('param') == len(parameter_names)
    assert line_2.count('qoi') == len(qoi_names)
    assert line_2.count('err') == len(error_names)
    assert line_2.count('qoi_v') == 0
    assert line_2.count('err_v') == 0

def test__write_header_section():

    cleanup_test()

    parameter_names = ['param{}'.format(i+1) for i in range(3)]
    qoi_names = ['qoi{}'.format(i+1) for i in range(5)]
    error_names = ['err{}'.format(i+1) for i in range(5)]
    
    datafile = PyposmatDataFile()
    datafile.write_header_section(
            parameter_names = parameter_names,
            qoi_names = qoi_names,
            error_names = error_names,
            filename =  datafile_out_fn)

    assert os.path.isfile(datafile_out_fn)

    datafile_read = PyposmatDataFile()
    datafile_read.read(filename = datafile_out_fn)

    assert len(datafile_read.parameter_names) == len(parameter_names)
    for i,v in enumerate(parameter_names):
        assert datafile_read.parameter_names[i] == v

    assert len(datafile_read.qoi_names) == len(qoi_names)
    for i,v in enumerate(qoi_names):
        assert datafile_read.qoi_names[i] == v

    assert len(datafile_read.error_names) == len(qoi_names)
    for i,v in enumerate(error_names):
        assert datafile_read.error_names[i] == v

    cleanup_test()

def test__write_simulation_results__no_filename():

    cleanup_test()

    parameter_names = ['param{}'.format(i+1) for i in range(3)]
    qoi_names = ['qoi{}'.format(i+1) for i in range(5)]
    error_names = ['err{}'.format(i+1) for i in range(5)]
    
    datafile = PyposmatDataFile()
    datafile.write_header_section(
            parameter_names = parameter_names,
            qoi_names = qoi_names,
            error_names = error_names,
            filename =  datafile_out_fn)

    sim_id = "test_id"
    
    results = OrderedDict()
    results['parameters'] = OrderedDict([(v,1.) for v in parameter_names])
    results['qois'] = OrderedDict([(v,2.) for v in qoi_names])
    results['errors'] = OrderedDict([(v,3.0) for v in error_names])
    
    datafile.write_simulation_results(
            sim_id,
            results)

    assert os.path.isfile(datafile_out_fn)

    datafile_read = PyposmatDataFile()
    datafile_read.read(filename=datafile_out_fn)

    

if __name__ == "__main__":

    parameter_names = ['param{}'.format(i+1) for i in range(3)]
    qoi_names = ['qoi{}'.format(i+1) for i in range(5)]
    error_names = ['err{}'.format(i+1) for i in range(5)]
    
    datafile = PyposmatDataFile()
    s = datafile.get_header_string(
            w_cluster_id = False,
            parameter_names = parameter_names,
            qoi_names = qoi_names,
            error_names = error_names)
    print(80*'=')
    print("{:^80}".format("header string"))
    print(80*'=')
    print(s)
    print(80*'=')
    print("{:^80}".format("attributes"))
    print(80*'=')
    print('sim_id:{}'.format(datafile.parameter_names))


