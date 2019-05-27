import os,copy,argparse
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataAnalyzer

if __name__ == "__main__":
    _fn_config=os.path.join("resources","pyposmat.config.in")
    _fn_data=os.path.join("resources","pyposmat.results.0.out")
    _fn_pareto_out=os.path.join("pyposmat.pareto.out")

    pda = PyposmatDataAnalyzer(fn_config=_fn_config,fn_data=_fn_data)
    pareto_df = pda.calculate_pareto_set()
    
    datafile=PyposmatDataFile()
    datafile.df = pareto_df
    datafile.parameter_names = pda.parameter_names
    datafile.qoi_names = pda.qoi_names
    datafile.error_names = pda.error_names
    datafile.names = ['sim_id'] \
            +datafile.parameter_names\
            +datafile.qoi_names\
            +datafile.error_names
    datafile.types = ['sim_id']\
            +len(datafile.parameter_names)*['param']\
            +len(datafile.qoi_names)*['qoi_names']\
            +len(datafile.error_names)*['error_names']
    datafile.write(_fn_pareto_out)
    datafile.read(_fn_pareto_out)
