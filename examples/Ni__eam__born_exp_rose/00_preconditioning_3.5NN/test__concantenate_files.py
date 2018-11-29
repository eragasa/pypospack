import pytest

import os
from collections import OrderedDict
from pathlib import Path
import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataFile

pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
print(pypospack_root_dir)

N_ranks = 128

sim_id_fmt = "{}_{}"
i_iteration = 1
i_sim_id = 0

pyposmat_data_out_fn = "data/pyposmat.results.{}.out".format(i_iteration)
pyposmat_data_filenames = ["rank_{}/pyposmat.results.out".format(i) for i in range(N_ranks)]
pyposmat_data_out = None
for i,v in enumerate(pyposmat_data_filenames):
    data_incremental = PyposmatDataFile()
    data_incremental.read(v)
    for row in data_incremental.df.iterrows():
        if pyposmat_data_out is None:
            pyposmat_data_out = PyposmatDataFile()
            pyposmat_data_out.write_header_section(
                    parameter_names=data_incremental.parameter_names,
                    qoi_names=data_incremental.qoi_names,
                    error_names=data_incremental.error_names,
                    filename=pyposmat_data_out_fn)
        sim_id = sim_id_fmt.format(i_iteration,i_sim_id)
        
        results=OrderedDict()
        results['parameters'] = OrderedDict()
        for p in pyposmat_data_out.parameter_names:
            results['parameters'][p] = row[1][p]
        
        results['qois'] = OrderedDict()
        for q in pyposmat_data_out.qoi_names:
            results['qois'][q] = row[1][q]
        
        results['errors'] = OrderedDict()
        for e in pyposmat_data_out.error_names:
            results['errors'][e] = row[1][e]

        pyposmat_data_out.write_simulation_results(
                sim_id=sim_id,
                results=results,
                cluster_id=None,
                filename=pyposmat_data_out_fn)
        i_sim_id += 1
