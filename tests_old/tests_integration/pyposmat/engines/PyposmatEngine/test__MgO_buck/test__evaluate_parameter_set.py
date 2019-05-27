import pytest

import os
from collections import OrderedDict
import pypospack.utils
from pypospack.pyposmat.engines import PyposmatEngine
from pypospack.pyposmat.data import PyposmatDataFile

pypospack_root_directory = pypospack.utils.get_pypospack_root_directory()
config_fn = os.path.join(pypospack_root_directory,'data/MgO_pareto_data/pyposmat.config.in')
simulation_directories = [
        'MgO_NaCl.lmps_min_all','MgO_NaCl.lmps_elastic','MgO_NaCl_001s.lmps_min_pos',
        'MgO_NaCl_fr_a.lmps_min_pos','MgO_NaCl_fr_c.lmps_min_pos','MgO_NaCl_sch.lmps_min_pos']

from pypospack.potentials.MgO import MgO_LewisCatlow
parameters=MgO_LewisCatlow['parameters']

def cleanup_simulation_directories():
    for d in simulation_directories:
        try:
           shutil.rmtree(d)
        except FileNotFoundError as e:
            pass

def test__evaluate_parameter_set():

    engine = PyposmatEngine(filename_in=config_fn,filename_out='pyposmat.config.out')
    engine.configure()
    results = engine.evaluate_parameter_set(parameters=parameters)
   
    assert 'parameters' in results
    assert 'qois' in results
    assert type(results['parameters']) is OrderedDict
    assert type(results['qois']) is OrderedDict
    assert all([v is not None for k,v in results['parameters'].items()])
    assert all([v is not None for k,v in results['qois'].items()])

def test__evaluate_tasks_through_the_task_mananager():
    engine = PyposmatEngine(filename_in=config_fn,filename_out='pyposmat.config.out')
    engine.configure()

    assert type(parameters) is OrderedDict
    assert type(engine.configuration.potential) is OrderedDict

    engine.task_manager.evaluate_tasks(parameters=parameters,potential=engine.configuration.potential)

    assert type(engine.task_manager.results) is OrderedDict

def test__send_results_to_the_qoi_manager():
    engine = PyposmatEngine(filename_in=config_fn,filename_out='pyposmat.config.out')
    engine.configure()
    engine.task_manager.evaluate_tasks(parameters=parameters,potential=engine.configuration.potential)
    engine.qoi_manager.calculate_qois(task_results=engine.task_manager.results)


if __name__ == "__main__":
    engine = PyposmatEngine(filename_in=config_fn,filename_out='pyposmat.config.out')
    engine.configure()
    results = engine.evaluate_parameter_set(parameters=parameters)
    print(results)
    if False:
        #print_parameters(parameters=parameters)
        m = [80*'=','{:^80}'.format('Potential Description'),80*'=']
        m += ['potential_description']
        for k,v in engine.configuration.potential.items():
            m += ['\t{} = {}'.format(k,v)]
        m += ['parameters']
        m += ['\t{} = {}'.format(k,v) for k,v in parameters.items()]
        m += [80*'=','{:80}'.format('Evaluate Parameter Set'),80*'=']
        engine.task_manager.evaluate_tasks(parameters=parameters,potential=engine.configuration.potential)
        m += ['engine.task_manager.results:']
        m += ['\ttype({})'.format(str(type(engine.task_manager.results)))]
        m += ['\t{} = {}'.format(k,v) for k,v in engine.task_manager.results.items()]
        for n_qoi, o_qoi in engine.qoi_manager.obj_Qoi.items():
            print(n_qoi)
            o_qoi.calculate_qois(task_results = engine.task_manager.results)
        print("\n".join(m))
        engine.qoi_manager.calculate_qois(task_results=engine.task_manager.results)


    
