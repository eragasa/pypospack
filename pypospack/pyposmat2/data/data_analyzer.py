# -*- coding: utf-8 -*-
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2016,2017"
__license__ = "Simplified BSD License"
__version__ = "1.0"

import os, copy
import numpy as np
import pandas as pd
from collections import OrderedDict
from pypospack.pyposmat2.data import PyposmatConfigurationFile
from pypospack.pyposmat2.data import PyposmatDataFile
import pypospack.pareto as pareto

class PyposmatDataAnalyzer(object):
    def __init__(self):
        self._pyposmat_configuration = None

    @property
    def pyposmat_configuration(self):
        return self._pyposmat_configuration

    @pyposmat_configuration.setter
    def pypospack_configuration(self,configuration):
        assert type(configuration) is PyposmatConfigurationFile

    @property
    def parameter_names(self): 
        return self._pyposmat_datafile.parameter_names

    @property
    def qoi_names(self):
        return self._pyposmat_datafile.qoi_names

    @property
    def error_names(self):
        return self._pyposmat_datafile.error_names

    @property
    def df(self):
        return self._df

    def read_configuration_file(self,filename):
        self._pyposmat_configuration = PyposmatConfigurationFile()
        self._pyposmat_configuration.read(filename)

    def read_data_file(self,filename):
        self._pyposmat_datafile = PyposmatDataFile(filename)
        self._pyposmat_datafile.read()

        self._df = copy.deepcopy(self._pyposmat_datafile.df)

    def calculate_pareto_set(self,df=None):
        if df is None:
            _df = copy.deepcopy(self.df)
        else:
            _df = copy.deepcopy(df)

        _df[self.error_names] = _df[self.error_names].abs()
        _results = []
        for i_row, row in _df.iterrows():
            _results.append([
                    i_row,
                    row[self.parameter_names].values.tolist(),
                    row[self.error_names].values.tolist()
                    ])
        is_pareto_idx = pareto.pareto(_results)
        print('n_pareto:{}'.format(len(is_pareto_idx)))
        if df is None:
            self._df['is_pareto'] = 0
            self._df.loc[is_pareto_idx,'is_pareto'] = 1
            return self._df
        else:
            _df = copy.deepcopy(df)
            _df['is_pareto'] = 0
            _df.loc[list(is_pareto_idx),'is_pareto'] = 1
            return _df
    
    def filter_performance_requirements(self,df=None):
        if df is None:
            _df = copy.deepcopy(self.df)
        else:
            _df = copy.deepcopy(df)
        
        _df[self.error_names] = _df[self.error_names].abs()
       
        _qoi_constraints = copy.deepcopy(
                self.pyposmat_configuration.qoi_constraints
            )

        is_survive_idx = []
        for k,v in _qoi_constraints.items():
            if k in self.error_names:
                is_survive_idx.append(_df.index[_df[k] < v])
                print('{} < {}: {}'.format(
                        k,
                        v,
                        len(_df.index[_df[k] < v])
                    ))
        is_survive_sets = [set(v) for v in is_survive_idx]
        is_survive_idx = list(set.intersection(*is_survive_sets))
        print('n_survive:{}'.format(len(is_survive_idx)))
        if df is None:
            self._df['is_survive'] = 0
            if len(is_survive_idx) > 0: 
                self._df.loc[is_survive_idx,'is_survive'] = 1
            return self._df
        else:
            _df = copy.deepcopy(df)
            _df.reset_index(drop=True)
            _df['is_survive'] = 0
            if len(is_survive_idx) > 0: 
                _df.loc[is_survive_idx,'is_survive'] = 1
            return _df
    def write_kde_file(self,filename):
        
        # <------ BEGIN OLD WORKING CODE
        # kde_df = self._df.loc[(self._df['is_pareto'] == 1) & (self._df['is_survive'] == 1)]
        # <------ END OF OLD WORKING CODE

        # filter by performance requirements first
        kde_df = None
        if len(self.pyposmat_configuration.qoi_constraints) > 0:
            kde_df = self.filter_performance_requirements(self._df)
            kde_df = kde_df.loc[kde_df['is_survive'] == 1]
            kde_df = kde_df.reset_index(drop=True)
        else:
            kde_df = copy.deepcopy(self._df)
    
        kde_df = self.calculate_pareto_set(df=kde_df)
        (n_rows_kde,n_cols_kde) = kde_df.shape
        print('n_samples_in_kde:{}'.format(n_rows_kde))
        print(kde_df)
        names = ['sim_id'] + self.parameter_names + self.qoi_names\
                + self.error_names
        types = ['sim_id'] \
                + len(self.parameter_names)*['param']\
                + len(self.qoi_names)*['qoi']\
                + len(self.error_names)*['err']

        str_list = []
        str_list.append(','.join(names))
        str_list.append(','.join(types))
        for i_row, row in kde_df[names].iterrows():
            str_list.append(','.join([str(v) for v in row.values.tolist()]))

        with open(filename,'w') as f:
            f.write("\n".join(str_list))
if __name__ == "__main__":
    
    data_directory = 'data'
    pyposmat_data_filename = 'pypospack.results.out'
    pyposmat_configuration_filename = 'pypospack.config.in'  
    data_analyzer = PypospackDataAnalyzer()
    data_analyzer.read_configuration_file(
            filename=pyposmat_configuration_filename)
    data_analyzer.read_data_file(
            filename=os.path.join(
                data_directory,
                pyposmat_data_filename))
    data_analyzer.calculate_pareto_set()
    data_analyzer.write_kde_file(
            filename=os.path.join(
                data_directory,
                'pyposmat.kde.out'))

    exit()
    datafile = PyposmatDataFile(filename=os.path.join(
        data_directory,pyposmat_data_filename))

    datafile.read()

    print(datafile.parameter_names)
    print(datafile.qoi_names)
    print(datafile.error_names)
    #print(datafile.df)

    df = copy.deepcopy(datafile.df)

    p_names = datafile.parameter_names
    q_names = datafile.qoi_names
    e_names = datafile.error_names

    p_column_idx = [df.columns.get_loc(n) for n in p_names if n in df.columns]
    q_column_idx = [df.columns.get_loc(n) for n in q_names if n in df.columns]
    e_column_idx = [df.columns.get_loc(n) for n in e_names if n in df.columns]

    df['sim_id'] = df['sim_id'].apply(np.int64)
    df[e_names] = df[e_names].abs()

    values = []
    i_iter = 0
    for i_row, row in df.iterrows():
        values.append(
                [
                    "{}_{}".format(i_iter,int(row['sim_id'])),
                    row[p_column_idx].values.tolist(),
                    row[e_column_idx].values.tolist()
                ])
    for idx,p,v in values:
        print('idx:',idx)
        print('p:',p)
        print('v:',v)
    
    is_pareto_idx = pareto.pareto([v for idx,p,v in values])
    df = copy.deepcopy(datafile.df)
    df['is_pareto'] = 0
    df.loc[is_pareto_idx,'is_pareto'] = 1
    n_results = len(values)
    n_pareto = len(is_pareto_idx)
    print('n_results={}'.format(n_results))
    print('n_pareto={}'.format(n_pareto))

    kde_parameters = df[df['is_pareto'] == 1][p_names]
    print(kde_parameters)
