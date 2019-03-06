import pytest

import os
from collections import OrderedDict

import pypospack.utils
from pypospack.pyposmat.data import PyposmatConfigurationFile

pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
configuration_filename = os.path.join(pypospack_root_dir,'data/MgO_pareto_data/pyposmat.config.in')
MgO = OrderedDict()
MgO['config_fn'] = os.path.join(pypospack_root_dir,'data/MgO_pareto_data/pyposmat.config.in')
MgO['n_iterations'] = 10 
def test__read():
    o = PyposmatConfigurationFile()
    o.read(filename=configuration_filename)

    assert type(o.n_iterations) is int
    assert o.n_iterations == 10

    assert type(o.qois) is OrderedDict
    for qoi_name,qoi_dict in o.qois.items():
        assert type(qoi_name) is str
        assert type(qoi_dict) is OrderedDict

        assert set(['qoi_type','structures','target']) == set(qoi_dict.keys())
        assert type(qoi_dict['qoi_type']) is str
        assert type(qoi_dict['structures']) is OrderedDict
        assert type(qoi_dict['target']) is float

    assert all([type(k) is str for k in o.qois])
    
def test__read__bad_pathname():
   o = PyposmatConfigurationFile()

   with pytest.raises(FileNotFoundError) as e:
       o.read(filename='badfilename')

def dev__read():
    o = PyposmatConfigurationFile()
    o.read(filename=configuration_filename)

    print('n_iterations({})={}'.format(type(o.n_iterations),o.n_iterations))
    print('qois({})'.format(type(o.qois)))
    for qoi_name,qoi_dict in o.qois.items():
        print('{}'.format(qoi_name))
        print('{}'.format(qoi_dict))
        
        
if __name__ == "__main__":
    dev__read()
