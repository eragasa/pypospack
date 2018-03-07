# -*- coding: utf-8 -*-
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2016,2017"
__license__ = "Simplified BSD License"
__version__ = "1.0"

import os, copy
import numpy as np
import pandas as pd
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile
import pypospack.pareto as pareto

class PyposmatDataAnalyzer(object):
    def __init__(self):
        self._configuration = None
        self._datafile = None
    @property
    def configuration(self):
        return self._configuration

    @configuration.setter
    def configuration(self,configuration):
        assert type(configuration) is PyposmatConfigurationFile

    @property
    def parameter_names(self):
        return self._datafile.parameter_names

    @property
    def qoi_names(self):
        return self._datafile.qoi_names

    @property
    def error_names(self):
        return self._datafile.error_names

    @property
    def df(self):
        return self._df

    def read_configuration_file(self,filename):
        self._configuration = PyposmatConfigurationFile()
        self._configuration.read(filename)

    def read_data_file(self,filename):
        self._datafile = PyposmatDataFile(filename)
        self._datafile.read()

        self._df = copy.deepcopy(self._datafile.df)

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
                self.configuration.qoi_constraints['filter_by_qoi_error']
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
    
    def calculate_d_metric(self,df):
        _df = copy.deepcopy(df)

        _qois = self.configuration.qois
        _error_names = self.error_names
        for i,en in enumerate(self.error_names):
            qn = self.qoi_names[i]
            q_target = _qois[qn]['target']
            _df[en] = np.abs(_df[en]/q_target)
            _df[en] = np.square(_df[en])

        _df['d_metric'] = np.sqrt(_df[_error_names].sum(axis=1))
        _df[_error_names] = df[_error_names]
        return _df

    def write_kde_file(self,filename):
        _qoi_constraints = self.configuration.qoi_constraints
        kde_df = copy.deepcopy(self._df)
        for k,v in _qoi_constraints.items():
            if k == 'filter_by_qoi_error':
                kde_df = self.filter_performance_requirements(kde_df)
                kde_df = kde_df.loc[kde_df['is_survive'] == 1]
                kde_df = kde_df.reset_index(drop=True)
            if k == 'filter_by_pareto':
                if v == True:
                    kde_df = self.calculate_pareto_set(df=kde_df)
                    kde_df = kde_df.loc[kde_df['is_pareto'] == 1]
                    kde_df = kde_df.reset_index(drop=True)
            if k == 'filter_by_dmetric':
                    (nr,nc) = kde_df.shape
                    kde_df = self.calculate_d_metric(kde_df)    
                    (nr,nc) = kde_df.shape
                    if v[1] == 'pct':
                        pct_to_keep = v[0]/100
                        n = int((nr * pct_to_keep) //1)
                    else:
                        n = min(v[0],nr)
                    kde_df = kde_df.nsmallest(n,'d_metric')
                    for en in self.error_names:
                        print(en,kde_df[en].max())
                    (nr,nc) = kde_df.shape
            (nr,nc) = kde_df.shape
            print('after {}: {} remainings'.format(k,nr))
        names = ['sim_id'] \
                + self.parameter_names \
                + self.qoi_names\
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
    configuration_filename = 'pypospack.config.in'
    data_analyzer = PypospackDataAnalyzer()
    data_analyzer.read_configuration_file(
            filename=configuration_filename)
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
