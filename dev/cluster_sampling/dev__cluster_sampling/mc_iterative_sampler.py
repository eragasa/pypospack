import os,shutil,sys
import numpy as np
from mpi4py import MPI
from pypospack.pyposmat.data import PyposmatConfigurationFile
from mc_sampler_iterate_w_cluster import PyposmatIterativeSampler
from mc_sampler_iterate_w_cluster import PyposmatClusterSampler
if __name__ == "__main__":
    pyposmat_data_directory = 'data'
    pyposmat_filename_in = os.path.join(
            pyposmat_data_directory,
            'pyposmat.config.in'
    )
    pyposmat_datafile_in = os.path.join(
            pyposmat_data_directory,
            "pyposmat.cluster.0.out"
    )
    #------------------------------------------------------------------------------
    # RUN PYPOSMAT 
    #------------------------------------------------------------------------------

    pyposmat_app = PyposmatIterativeSampler(
        configuration_filename = pyposmat_filename_in)
    pyposmat_app.data_directory = pyposmat_data_directory
    pyposmat_app.read_configuration_file()
    pyposmat_app.run_all()
    exit()


    o = PyposmatClusterSampler()
    try:
        o.read_configuration_file(
                filename=pyposmat_filename_in)
    except FileNotFoundError as e:
        print("attempted to configure with")
        print("    pyposmat_configuration_fn={}".format(
            pyposmat_filename_in)
        )
        raise e
    print("o.pyposmat_configuration_fn={}".format(o.configuration_fn))
    print("type(o.configuration)={}".format(type(o.configuration)))
    assert type(o.configuration) is PyposmatConfigurationFile
    print("structure_directory={}".format(
        o.structure_directory))
    print("parameter_names=:")
    for p in o.parameter_names:
        print("    {}".format(p))
    print("qoi_names:")
    for q in o.qoi_names:
        print("   {}".format(q))
    print("error_names:")
    for e in o.error_names:
        print("   {}".format(e))

    # read the pyposmat datafile
    o.configure_pyposmat_datafile_in(
        filename=pyposmat_datafile_in
    )
    print(o.data)
    print(o.data.df)
    print(o.data.df.columns)
    i_iteration = 0
    #for k,v in o.configuration.sampling_type.items():
    #    print(k,v)
    o.run_simulations(i_iteration=0)



