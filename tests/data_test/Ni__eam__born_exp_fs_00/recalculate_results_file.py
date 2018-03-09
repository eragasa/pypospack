import os,copy
import pandas as pd
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatDataAnalyzer
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

_fn_config = os.path.join("","data__Ni__eam__born_exp_fs_01","pyposmat.config.in")
_fn_results_in = os.path.join("data__Ni__eam__born_exp_fs_01","pyposmat.results.0.out")
_fn_results_out = os.path.join("data__Ni__eam__born_exp_fs_01","pyposmat.results.0a.out")

config = PyposmatConfigurationFile()
config.read(_fn_config)
qoi_targets = config.qoi_targets
print(config.qoi_targets)
print(list(config.qois))
data_in = PyposmatDataFile()
data_in.read(_fn_results_in)
data_out = PyposmatDataFile()
data_out.parameter_names = data_in.parameter_names
data_out.qoi_names = list(config.qois)
data_out.error_names = ['{}.err'.format(q) for q in data_out.qoi_names]
data_out.names = ["sim_id"]\
        +data_out.parameter_names\
        +data_out.qoi_names\
        +data_out.error_names
data_out.types = ["sim_id"]\
        +len(data_out.parameter_names)*['param']\
        +len(data_out.qoi_names)*['qoi']\
        +len(data_out.error_names)*['err']

def calculate_bulk_modulus(c11,c12,c44):
    return (c11+2*c12)/3

def calculate_shear_modulus(c11,c12,c44):
    return (c11-c12)/2

data_out_lists = []
for i,row in data_in.df.iterrows():

    in_row_results = row.to_dict(into=OrderedDict)
    out_row_results = OrderedDict()

    for k in (["sim_id"]+data_out.parameter_names):
        out_row_results[k] = in_row_results[k]
    for k in (data_out.qoi_names):
        try:
            out_row_results[k] = in_row_results[k]
        except KeyError as e:
            if k == 'Ni_fcc.B':
                c11 = in_row_results['Ni_fcc.c11']
                c12 = in_row_results['Ni_fcc.c12']
                c44 = in_row_results['Ni_fcc.c44']
                out_row_results[k] = calculate_bulk_modulus(c11,c12,c44)
            elif k == 'Ni_fcc.G':
                c11 = in_row_results['Ni_fcc.c11']
                c12 = in_row_results['Ni_fcc.c12']
                c44 = in_row_results['Ni_fcc.c44']
                out_row_results[k] = calculate_bulk_modulus(c11,c12,c44)
            else:
                raise
    for k in (data_out.qoi_names):
        out_row_results["{}.err".format(k)] = out_row_results[k] - qoi_targets[k]
    data_out_lists.append([out_row_results[k] for k in data_out.names])
data_out.df=pd.DataFrame(data_out_lists,columns=data_out.names)
data_out.write(_fn_results_out)
