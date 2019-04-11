import copy,os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.potential import EamPotential
from pypospack.eamtools import EamSetflFile


if __name__ == "__main__":
    data_directory = "data/new_data"
    
    config_fn = os.path.join(data_directory, "pyposmat.config.in")
    config = PyposmatConfigurationFile()
    config.read(filename=config_fn)
    
    data_fn = os.path.join(data_directory, "pyposmat.kde.3.out")
    datafile=PyposmatDataFile()
    datafile.read(filename=data_fn)

    for iqn,qn in enumerate(datafile.qoi_names):
        en = "{}.err".format(qn)
        nen = "{}.nerr".format(qn)
        q = config.qoi_targets[qn]
        datafile.df[nen] = datafile.df[qn]/q-1
    (nrows,ncols) = datafile.df.shape

    normederr_names = ['{}.nerr'.format(q) for q in datafile.qoi_names]

    datafile.df['d_metric'] = np.sqrt(np.square(datafile.df[normederr_names]).sum(axis=1))
    subselect_df = datafile.df.nsmallest(1, "d_metric")
    subselect_df.set_index("sim_id")

    symbols = ["Ni"]
    func_pair_name = "bornmayer"
    func_density_name = "eam_dens_exp"
    func_embedding_name = "eam_embed_fs"

    eam_potential = EamPotential(symbols=symbols, func_pair=func_pair_name,
                                 func_density=func_density_name, func_embedding=func_embedding_name)
    
    parameters = OrderedDict()
    for name in datafile.parameter_names:
        parameters[name] = subselect_df[name].values[0]

    # plotting portion copied from pypospack.tests.tests_integration.eamtools.EamSetflFile.dev_EamSetflFile

    setfl_fn = "Ni__eam__born_exp_fs.eam.alloy"
    n_r = 2000
    r_max = 10.0
    r_cut = 10.0
    n_rho = 2000
    rho_max = 10.0

    eam_potential.write_setfl_file(filename=setfl_fn, symbols=symbols,
                                   Nr=n_r, rmax=r_max, rcut=r_cut,
                                   Nrho=n_rho, rhomax=rho_max, 
                                   parameters=parameters)

    setfl_file = EamSetflFile()
    setfl_file.read(setfl_fn)

    r = setfl_file.r
    rho = setfl_file.rho
    pair = setfl_file.func_pairpotential['Ni.Ni']
    embed = setfl_file.func_embedding['Ni']
    dens = setfl_file.func_density['Ni']

    fig, axes = plt.subplots(3,1)
    r_low = 1
    r_high = 10
    phi_low = pair.min()
    phi_high = 0.5
    axes[0].plot(r,pair/r)
    axes[0].set_xlim(r_low,r_high)
    axes[0].set_ylim(phi_low,phi_high)
    axes[1].plot(rho,embed)
    axes[2].plot(r,dens)
    axes[2].set_xlim(r_low,r_high)
    axes[2].set_ylim(0,10)
    plt.show()



