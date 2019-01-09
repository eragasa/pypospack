from pypospack.pyposmat.data import PyposmatDataAnalyzer

if __name__ == "__main__":
    import os
    import pypospack.utils

    config_fn = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'examples','Ni__eam__born_exp_rose','01_preconditioning_3.5NN','data',
            'pyposmat.config.in')

    data_fn = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'examples','Ni__eam__born_exp_rose','01_preconditioning_3.5NN','data',
            'pyposmat.results.0.out')

    # check if files exist
    if os.path.isfile(config_fn):
        print('config_fn:{}'.format(config_fn))
    else:
        print('cannot find config_fn:{}'.format(config_fn))
        exit()

    
    if os.path.isfile(data_fn):
        print('data_fn:{}'.format(data_fn))
    else:
        print('cannot find data_fn:{}'.format(data_fn))
        exit()

    # initialization
    o = PyposmatDataAnalyzer(fn_config=config_fn,fn_data=data_fn)
    print(o.parameter_names)
    results_nrows, results_ncols = o.data.df.shape
    # calculate_pareto_set
    pareto_df = o.calculate_pareto_set(df=o.data.df)
    pareto_nrows, pareto_ncols = pareto_df.loc[pareto_df["is_pareto"] == True].shape

    print("n_results:{}".format(results_nrows))
    print("n_pareto:{}".format(pareto_nrows))

    print(pareto_df["is_pareto"])
