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
from pypospack.pareto import pareto

def filter_by_pareto_set_membership(df,error_names):
    pass

def filter_by_znormalized_errors(df,percentile,qoi_names):
    for qn in qoi_names:
        en = "{}.err".format(qn)
        zen = "{}.zerr".format(qn)
        df[zen] = (df[en]-df[en].mean())/df[en].std()

    zerror_names = ["{}.zerr".format(q) for q in qoi_names]
    df['z_err_dist'] = np.sqrt(np.square(df[zerror_names]).sum(axis=1))

    nr0,nc0 = df.shape
    nr1 = int(percentile*nr0//1)
    df=df.nsmallest(nr1,"z_err_dist").reset_index(drop=True)

    return df

class PyposmatDataAnalyzer(object):
    def __init__(self,fn_config=None,fn_data=None):
        self._configuration = None
        self._data = None
        self.obj_log = None

        if fn_config is not None:
            self.fn_config = fn_config
            self.read_configuration_file(fn_config)
        if fn_data is not None:
            self.fn_data = fn_data
            self.read_data_file(fn_data)

    @property
    def configuration(self):
        return self._configuration

    @configuration.setter
    def configuration(self,configuration):
        assert type(configuration) is PyposmatConfigurationFile
        self._configuration = configuration

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self,data):
        self._data = data

    @property
    def datafile(self):
        return self._data

    @datafile.setter
    def datafile(self,datafile):
        assert type(datafile) is PyposmatDataFile
        self._data = datafile

    @property
    def parameter_names(self):
        return self._parameter_names

    @property
    def qoi_names(self):
        return self._qoi_names

    @property
    def error_names(self):
        return self._error_names

    @property
    def qoi_targets(self):
        return OrderedDict([(qn,qv['target'] )for qv,qv in self.configuration.qois.items()])

    @property
    def df(self):
        return self._df

    def log(self,m):
        if self.obj_log is None:
            print(m)

    def read_configuration(self,filename):
        self._configuration = PyposmatConfigurationFile()
        self._configuration.read(filename=filename)

        m = "configuration:{}".format(filename)
        self.log(m)

    def read_configuration_file(self,filename):
        self.read_configuration(filename=filename)

    def read_data(self,filename):
        self._data = PyposmatDataFile()
        self._data.read(filename=filename)

        m = "data:{}".format(filename)
        self.log(m)

        self._df = copy.deepcopy(self.data.df)
        self._parameter_names = list(self.data.parameter_names)
        self._error_names = list(self.data.error_names)
        self._qoi_names = list(self.data.qoi_names)
        self._score_names = None
    
    def read_data_file(self,filename):
        
        self.read_data(filename=filename)

    def calculate_pareto_set(self,df=None,v=None):
        
        if df is None:
            _df = self.data.df
        else:
            _df = copy.deepcopy(df)

        # define names for absolute errors
        abs_error_names = ["{}.abserr".format(k) for k in self.qoi_names]
        for i,qn in enumerate(self.qoi_names):
            en  = "{}.err".format(qn)
            aen = "{}.abserr".format(qn)
            _df[aen] = _df[en].abs()

        is_pareto_idx = pareto(_df[abs_error_names].values.tolist())

        _df['is_pareto'] = 0
        _df.loc[list(is_pareto_idx),'is_pareto'] = 1
        return _df

    def filter_performance_requirements(self,df=None):
        """
        Args:
            df (pandas.DataFrame)
        """
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
            q= self.qoi_targets[self.qoi_names[i]]

            _df[en] = np.abs(_df[en]/q)
            _df[en] = np.square(_df[en])

        _df['d_metric'] = np.sqrt(_df[_error_names].sum(axis=1))
        _df[_error_names] = df[_error_names]
        return _df

    def filter_with__qoi_constraints(self,kde_df,qoi_constraints):
        _df = copy.deepcopy(kde_df)
        m = [80*'-',"filtering by qoi_constraints",80*'-']
        self.log("\n".join(m))

        for qn in self.datafile.qoi_names:
            aen = "{}.abserr".format(qn)
            en = "{}.err".format(qn)
            _df[aen] = _df[en].abs()

        for qoic_n,qoic_v in qoi_constraints.items():
            nr0,nc0 =_df.shape
            if qoic_v[0] == '<':
                _df = _df[_df[qoic_n] < qoic_v[1]]
            elif qoic_v[0] == '>':
                _df = _df[_df[qoic_n] > qoic_v[1]]
            elif qoic_v[0] == '=':
                _df = _df[_df[qoic_n] == qoic_v[1]]
            else:
                raise ValueError('unknown operator, {}'.format(k))
            nr1,nc0 =_df.shape
            m = "{:>20} {:^3} {:<10} {:10} {:10} {:10}".format(qoic_n,qoic_v[0],qoic_v[1],nr1,nr0,nr0-nr1)
            self.log(m)
        return _df

    def filter_by_znormalized_errors(self,kde_df,qoi_constraints):
        percentile = qoi_constraints['percentile']
        return filter_by_znormalized_errors(
                df=kde_df,
                percentile=percentile,
                qoi_names=self.qoi_names)

    def write_kde_file(self,filename,qoi_constraints=None):
        """
        """

        if qoi_constraints is None:
            _qoi_constraints = self.configuration.qoi_constraints
        else:
            _qoi_constraints = qoi_constraints

        kde_df = copy.deepcopy(self._df)

        # filtering by the qoi constraints
        for k,v in _qoi_constraints.items():

            if k == 'qoi_constraints':
                kde_df = self.filter_with__qoi_constraints(kde_df,v)
            elif k == "filter_by__d_zerror":
                kde_df = self.filter_by_znormalized_errors(kde_df,v)
            elif k == 'filter_by_pareto' or k == 'select_pareto_only':
                kde_df = self.calculate_pareto_set(kde_df,v)
                kde_df = kde_df.loc[kde_df['is_pareto'] == 1]
            elif k == 'filter_by_dmetric':
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
            else:
                raise ValueError("unknown qoi_constraint method {}".format(k))
            
            kde_df = kde_df.reset_index(drop=True)
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
    pass
