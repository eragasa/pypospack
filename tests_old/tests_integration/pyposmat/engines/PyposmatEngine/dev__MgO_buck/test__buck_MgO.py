import pytest
import os
from collections import OrderedDict
import pypospack.utils
from pypospack.pyposmat.engines import PyposmatEngine
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

_pypospack_root_directory = pypospack.utils.get_pypospack_root_directory()
_pyposmat_config_fn = os.path.join(_pypospack_root_directory,'data','MgO_pareto_data','pyposmat.config.in')

def test__init():
    engine = PyposmatEngine(filename_in=_pyposmat_config_fn,filename_out='pyposmat.config.out')

def test__configure():
    engine = PyposmatEngine(filename_in=_pyposmat_config_fn,filename_out='pyposmat.config.out')
    engine.configure()

def test__create_base_directories():
    engine = PyposmatEngine(filename_in=_pyposmat_config_fn,filename_out='pyposmat.config.out')
    engine.create_base_directories()

def test__read_configuration_file():
    engine = PyposmatEngine(filename_in=_pyposmat_config_fn,filename_out='pyposmat.config.out')
    engine.read_configuration_file()

    assert isinstance(engine.configuration,PyposmatConfigurationFile)
    assert engine.configuration.qois is not None
def test__configure_qoi_manager():
    engine = PyposmatEngine(filename_in=_pyposmat_config_fn,filename_out='pyposmat.config.out')
    engine.read_configuration_file()
    engine.configure_qoi_manager()

def test__configure():
    # this test are the submethods in the configure() method
    engine = PyposmatEngine(filename_in=_pyposmat_config_fn,filename_out='pyposmat.config.out')
    engine.create_base_directories()
    engine.read_configuration_file()
    engine.configure_qoi_manager()
    engine.configure_task_manager()

if __name__ == "__main__":
    engine = PyposmatEngine(filename_in=_pyposmat_config_fn,filename_out='pyposmat.config.out')
    engine.read_configuration_file()
    #print(engine.configuration)
    #print(engine.configuration.read(filename=_pyposmat_config_fn))
    print(engine.configuration.configuration)

    if False:
        print("structure_directory:{}".format(engine.structure_directory))
        print("structures")
        for k,v in engine.structures['structures'].items():
            print(k,v)

        from pypospack.potentials.MgO import MgO_LewisCatlow
        _parameters=MgO_LewisCatlow['parameters']
        _results=engine.evaluate_parameter_set(parameters=_parameters)

        print(80*'-')
        print('parameters')
        for k,v in _results['parameters'].items():
            print(k,v)

        print(80*'-')
        print('qois')
        for k,v in _results['qois'].items():
            print(k,v)


