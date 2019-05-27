import os,copy,argparse
from pypospack.pyposmat.data import PyposmatDataAnalyzer

if __name__ == "__main__":
    _fn_config=os.path.join(
         "data__Ni__eam__born_exp_bjs_05",
         "pyposmat.config.in")
    _fn_data=os.path.join(
         "data__Ni__eam__born_exp_bjs_04",
         "pyposmat.kde.9.out")
    _fn_kde_out=os.path.join(
         "data__Ni__eam__born_exp_bjs_05",
         "pyposmat.kde.00.out")

    pda = PyposmatDataAnalyzer(
        fn_config=_fn_config,
        fn_data=_fn_data)
    pda.write_kde_file(filename=_fn_kde_out)
