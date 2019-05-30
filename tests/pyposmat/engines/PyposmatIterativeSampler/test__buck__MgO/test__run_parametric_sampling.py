import pytest 

import os,shutil,sys
import numpy as np
from mpi4py import MPI
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.engines import PyposmatIterativeSampler

pyposmat_data_dir = 'data'
config_fn = os.path.join(pyposmat_data_dir,'pyposmat.config.in')
i_iteration = 0

def test__run_parametric_sampling():
    o = PyposmatIterativeSampler(configuration_filename=config_fn)
    o.data_dir = pyposmat_data_dir
    o.read_configuration_file()

    # ensure that the configuration is actually setup as parametric
    assert o.configuration.sampling_type[i_iteration]['type'] == 'parametric'

    #this setup happens in run_all()
    o.setup_mpi_environment()
    o.initialize_data_directory()
    o.i_iteration = i_iteration
    o.log_iteration_information(i_iteration=o.i_iteration)
    
    #this setup happens in run_simulations()
    o.initialize_rank_directory()

    # ensure that the paths are absolute paths
    assert os.path.isabs(o.rank_directory)
    assert os.path.isabs(o.configuration_filename)
    config_filename = o.configuration_filename
    results_filename = os.path.join(o.rank_directory,'pyposmat.results.{}.out')
    bad_parameters_filename = os.path.join(o.rank_directory,'pyposmat.bad_parameters.{}.out')

    os.chdir(o.rank_directory)

    o.determine_rv_seeds()

    o.initialize_sampler(
            config_fn=config_filename,
            results_fn=results_filename,
            mpi_rank=o.mpi_rank,
            mpi_size=o.mpi_size,
            o_log=o.o_log)


    o.run_parametric_sampling(i_iteration=i_iteration)

def dev__run_parametric_sampling():
    o = PyposmatIterativeSampler(configuration_filename=config_fn)
    o.data_dir = pyposmat_data_dir
    o.read_configuration_file()

    # ensure that the configuration is actually setup as parametric
    assert o.configuration.sampling_type[i_iteration]['type'] == 'parametric'

    #this setup happens in run_all()
    o.setup_mpi_environment()
    o.initialize_data_directory()
    o.i_iteration = i_iteration
    o.log_iteration_information(i_iteration=o.i_iteration)
    
    #this setup happens in run_simulations()
    o.initialize_rank_directory()

    # ensure that the paths are absolute paths
    assert os.path.isabs(o.rank_directory)
    assert os.path.isabs(o.configuration_filename)
    config_filename = o.configuration_filename
    results_filename = os.path.join(o.rank_directory,'pyposmat.results.{}.out')
    bad_parameters_filename = os.path.join(o.rank_directory,'pyposmat.bad_parameters.{}.out')

    os.chdir(o.rank_directory)

    o.determine_rv_seeds()

    o.initialize_sampler(
            config_fn=config_filename,
            results_fn=results_filename,
            mpi_rank=o.mpi_rank,
            mpi_size=o.mpi_size,
            o_log=o.o_log)

    o.run_parametric_sampling(i_iteration=i_iteration)
if __name__ == "__main__":
    dev__run_parametric_sampling()
