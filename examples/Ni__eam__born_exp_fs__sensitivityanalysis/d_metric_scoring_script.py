import pandas as pd
import numpy as np
from pypospack.pyposmat.data import PyposmatDataFile, PyposmatDataAnalyzer


if __name__ == "__main__":

    data_file_name = "/Users/seaton/python-repos/pypospack/examples/Ni__eam__born_exp_fs__sensitivityanalysis/data__from_pareto_optiization/pyposmat.results.5.out"
    config_file_name = "/Users/seaton/python-repos/pypospack/examples/Ni__eam__born_exp_fs__sensitivityanalysis/data__from_pareto_optiization/pyposmat.config.in"
    data = PyposmatDataFile(filename=data_file_name)
    data.read(data_file_name)

    analyzer = PyposmatDataAnalyzer()
    analyzer.read_configuration_file(filename=config_file_name)
    analyzer.read_data_file(filename=data_file_name)

    df = analyzer.calculate_d_metric(df=analyzer.datafile.df)

    data.df = df
    data.subselect_by_score(score_name="d_metric", n=100)
    # print(data.sub_df)

    param_stdev_df = data.sub_parameter_df.std(axis=0)
    param_mean_df = data.sub_parameter_df.mean(axis=0)
    print("parameter standard deviations:\n{}".format(param_stdev_df))
    print("parameter means:\n{}".format(param_mean_df))

    data.subselect_by_score(score_name="d_metric", n=1)
    print("best parameterization by d_metric:\n{}".format(data.sub_parameter_df))
