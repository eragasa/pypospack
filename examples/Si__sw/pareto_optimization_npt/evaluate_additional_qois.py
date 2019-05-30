import os,shutil,sys,copy
import numpy as np
import pandas as pd
from mpi4py import MPI
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.engines import PyposmatFileSampler

def check_reference_potentials(o_sampler):
    print(80*'-')
    print('check reference potentials')
    print(80*'-')
    reference_potential_names = ["LC","BG1","BG2"]
    for v in reference_potential_names:
        print("{}:{}".format(
            v,
            v in o_sampler.reference_potentials
            ))

def print_qoi_database(o_sampler):
    print(80*'-')
    print("{:^80}".format('QOI DATABASE FOR THIS ITERATION'))
    print(80*'-')
    for qoi_name,qoi_config in o_sampler.qoi_manager.qoidb.qois.items():
        print(qoi_name)
        for k,v in qoi_config.items():
            print("\t{}:{}".format(k,v))

def print_qoi_manager(o_sampler):
    print(80*'-')
    print("{:^80}".format('LIST OF TASKS DETERMINED BY QOIMANAGER'))
    print(80*'-')

    for task_name,task_configuration in o_sampler.qoi_manager.tasks.items():
        print("{}".format(task_name))
        for k,v in task_configuration.items():
            print("\t{}:{}".format(k,v))

def print_task_manager(o_sampler):

    print("\to_sampler.task_manager:{}".format(type(o_sampler.task_manager)))
    print("\to_sampler.task_manager.obj_Task:{}".format(type(o_sampler.task_manager.obj_Task)))
    for k,v in o_sampler.task_manager.obj_Task.items():
        print("\t\t{}:{}".format(k,type(v)))

def setup__configure_qoi_manager():
    config_directory = "./data"
    config_fn = os.path.join(config_directory,'pyposmat.config.in')
    
    data_directory = "../../data/MgO_pareto_data"
    datafile_in_fn = os.path.join(data_directory,'culled_005.out')

    output_directory = "./"
    datafile_out_fn = os.path.join(output_directory,'qoiplus_005.out')
    
    o_sampler = PyposmatFileSampler(
            config_fn=config_fn,
            data_in_fn=datafile_in_fn,
            data_out_fn=datafile_out_fn)
    o_sampler.create_base_directories()
    o_sampler.read_configuration_file()

def test__configure_qoi_manager__from_OrderedDict():
    setup__configure_qoi_manager()

    _qois = OrderedDict()
    for k,v in o_sampler.configuration.qois_validation.items():
        _qois[k] =v

    o_sampler.configure_qoi_manager(qois=_qois)

    from pypospack.qoi import QoiManager,QoiDatabase
    assert type(o_sampler.qoi_manager) is QoiManager
    assert type(o_sampler.qoi_manager.qoidb) is QoiDatabase

    for k,v in _qois.items():
        assert k in o_sampler.qoi_manager.qoi_db
    
def test__configure_task_manager():
    o_sampler = PyposmatFileSampler(
            config_fn=config_fn,
            data_in_fn=datafile_in_fn,
            data_out_fn=datafile_out_fn)
    o_sampler.create_base_directories()
    o_sampler.read_configuration_file()
    o_sampler.configure_qoi_manager(use_fitting_qois=False,use_testing_qois=True)
    o_sampler.configure_task_manager()
    
    from pypospack.task import TaskManager
    assert type(o_sampler.task_manager) is TaskManager

def test__subselect_by_dmetric():
    
    o_sampler = PyposmatFileSampler(
            config_fn=config_fn,
            data_in_fn=datafile_in_fn,
            data_out_fn=datafile_out_fn
            )
    o_sampler.create_base_directories()
    o_sampler.read_configuration_file()
    o_sampler.configure_qoi_manager(use_fitting_qois=False,use_testing_qois=True)
    o_sampler.configure_task_manager()
    o_sampler.configure_datafile_out()
    o_sampler.subselect_by_dmetric(nsmallest=n_smallest)
    
    import pandas as pd
    assert type(o_sampler.subselect_df) is pd.DataFrame
if __name__ == "__main__":

    n_smallest = 100

    config_directory = "./data"
    config_fn = os.path.join(config_directory,'pyposmat.config.in')
    
    data_directory = "../../data/MgO_pareto_data"
    datafile_in_fn = os.path.join(data_directory,'culled_005.out')

    output_directory = "./"
    datafile_out_fn = os.path.join(output_directory,'qoiplus_005.out')

    from pypospack.pyposmat.data import PyposmatConfigurationFile
    o_config=PyposmatConfigurationFile()
    o_config.read(filename=config_fn)

    o_sampler = PyposmatFileSampler(
            config_fn=config_fn,
            data_in_fn=datafile_in_fn,
            data_out_fn=datafile_out_fn
            )

    check_reference_potentials(o_sampler=o_sampler)

    o_sampler.create_base_directories()
    o_sampler.read_configuration_file()
    
    # Determine which QOIS you want to calculate
    # calculate only the fitting qois
    #o_sampler.configure_qoi_manager(use_fitting_qois=True,use_testing_qois=False)
    # Calculate only the testing qois
    o_sampler.configure_qoi_manager(use_fitting_qois=False,use_testing_qois=True)
    # Calculate all qois
    #o_sampler.configure_qoi_manager(use_fitting_qois=True,use_testing_qois=True)

    print_qoi_database(o_sampler=o_sampler)
    print_qoi_manager(o_sampler=o_sampler)
    
    o_sampler.configure_task_manager()
    print_task_manager(o_sampler=o_sampler)
    

    o_sampler.configure_datafile_out()
    o_sampler.subselect_by_dmetric(nsmallest=n_smallest)
   
    # these have been tested to work
    #o_sampler._sample_from_subselect_df()
    #o_sampler._sample_from_reference_potentials() 
   
    o_sampler.run_simulations()
    
    exit()
