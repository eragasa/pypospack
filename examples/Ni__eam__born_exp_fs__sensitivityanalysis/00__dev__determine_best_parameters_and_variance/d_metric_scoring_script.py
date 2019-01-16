import pandas as pd
import numpy as np
from pypospack.pyposmat.data import PyposmatDataFile, PyposmatDataAnalyzer


if __name__ == "__main__":
    import os
    import pypospack.utils
    _pypospack_root = pypospack.utils.get_pypospack_root_directory()
    _data_in_directory = os.path.join(_pypospack_root,'examples','Ni__eam__born_exp_fs__sensitivityanalysis','data__from_pareto_optimization')
    _pyposmat_data_fn = os.path.join(_data_in_directory,'pyposmat.kde.6.out')
    _pyposmat_config_fn = os.path.join(_data_in_directory,'pyposmat.config.in')


    analyzer = PyposmatDataAnalyzer()
    analyzer.read_configuration_file(filename=_pyposmat_config_fn)
    analyzer.read_data_file(filename=_pyposmat_data_fn)
    df = analyzer.calculate_d_metric(df=analyzer.datafile.df)

    data = PyposmatDataFile()
    data.read(filename=_pyposmat_data_fn)
    data.df = df
    data.subselect_by_score(score_name="d_metric", n=100)
    # print(data.sub_df)

    param_stdev_df = data.sub_parameter_df.std(axis=0)
    param_mean_df = data.sub_parameter_df.mean(axis=0)
    print("parameter standard deviations:\n{}".format(param_stdev_df))
    print("parameter means:\n{}".format(param_mean_df))

    data.subselect_by_score(score_name="d_metric", n=1)
    print("best parameterization by d_metric:\n{}".format(data.sub_parameter_df))
