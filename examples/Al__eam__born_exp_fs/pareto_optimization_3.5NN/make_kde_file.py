import os
from pypospack.pyposmat.data import PyposmatDataFile, PyposmatConfigurationFile
import numpy as np


if __name__ == "__main__":
    data_fn = "../preconditioning_3.5NN/data/pyposmat.kde.3.out"
    config_fn = "../preconditioning_3.5NN/data/pyposmat.config.in"

    o_data = PyposmatDataFile()
    o_data.read(filename=data_fn)
    o_config = PyposmatConfigurationFile()
    o_config.read(filename=config_fn)

    a0_max = o_config.qoi_targets["Al_fcc.a0"] + 0.5
    a0_min = o_config.qoi_targets["Al_fcc.a0"] - 0.5
    e_coh_max = o_config.qoi_targets["Al_fcc.E_coh"] + 2
    e_coh_min = o_config.qoi_targets["Al_fcc.E_coh"] - 2

    print("n points initial: {}".format(len(o_data.df)))
    o_data.df = o_data.df[o_data.df["Al_fcc.a0"] > a0_min]
    o_data.df = o_data.df[o_data.df["Al_fcc.a0"] < a0_max]
    o_data.df = o_data.df[o_data.df["Al_fcc.E_coh"] > e_coh_min]
    o_data.df = o_data.df[o_data.df["Al_fcc.E_coh"] < e_coh_max]
    print("n points final: {}".format(len(o_data.df)))

    print()
    print("a0 target: {}".format(o_config.qoi_targets["Al_fcc.a0"]))
    print("a0 min: {}".format(o_data.df["Al_fcc.a0"].min()))
    print("a0 max: {}".format(o_data.df["Al_fcc.a0"].max()))
    print()
    print("E_coh target: {}".format(o_config.qoi_targets["Al_fcc.E_coh"]))
    print("E_coh min: {}".format(o_data.df["Al_fcc.E_coh"].min()))
    print("E_coh max: {}".format(o_data.df["Al_fcc.E_coh"].max()))

    out_fn = "./data/pyposmat.kde.0.out"
    o_data.write(filename=out_fn)
