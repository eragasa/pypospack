import os,copy,argparse
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataAnalyzer

def str__qoi_constraints(qoi_constraints):
    s = []

    for k,v in qoi_constraints.items():
        if k == 'qoi_constraints':
            s.append(k)
            for qoic_n,qoic_v in v.items():
                s.append("{:5} {:>20} {:^3} {:<10}".format("",qoic_n,qoic_v[0],qoic_v[1]))
        else:
            s.append(",".join([str(k),str(v)]))
    return "\n".join(s)

def print__qoi_constraints(qoi_constraints):
    s = str__qoi_constraints(qoi_constraints)
    print(s)

class Dev__PyposmatDataAnalyzer(PyposmatDataAnalyzer):
    pass

if __name__ == "__main__":
    _fn_config=os.path.join("data__Ni__eam__born_exp_fs_01","pyposmat.config.in")
    _fn_data=os.path.join("data__Ni__eam__born_exp_fs_01","pyposmat.results.0a.out")
    _fn_kde_out=os.path.join("pyposmat.kde.out")

    pda = PyposmatDataAnalyzer(
        fn_config=_fn_config,
        fn_data=_fn_data)
    pda.write_kde_file(filename=_fn_kde_out)
