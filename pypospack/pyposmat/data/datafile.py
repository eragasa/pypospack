import os,copy
from collections import OrderedDict
import pandas as pd
import numpy as np

class PyposmatDataFile(object):

    def __init__(self,filename=None):
        self.filename = filename

        self.names = None
        self.types = None
        self.parameter_names = None
        self.qoi_names = None
        self.error_names = None
        self.score_names = None
        self.qoi_references = None
        self.scaling_factors = None

        self.df = None
        self.parameter_df = None
        self.error_df = None
        self.qoi_df = None
        self.rescaled_error_df = None

        self.sub_indices = None
        self.sub_df = None
        self.sub_parameter_df = None
        self.sub_qoi_df = None
        self.sub_error_df = None

    def write_header_section(self,
            parameter_names,
            qoi_names,
            error_names,
            filename=None):

        if filename is not None:
            assert isinstance(filename,str)
            self.filename=filename

        assert isinstance(parameter_names,list)
        assert isinstance(qoi_names,list)
        assert isinstance(error_names,list)

        #_a0 = self.configuration['potential']['a0']
        self.parameter_names = list(parameter_names)
        self.qoi_names = list(qoi_names)
        self.error_names = list(error_names)

        self.names = ['sim_id']\
                + list(parameter_names)\
                + list(qoi_names)\
                + list(error_names)
        self.types = ['sim_id']\
                + len(parameter_names) * ['param']\
                + len(qoi_names) * ['qoi']\
                + len(error_names) * ['err']\

        _header_str = ",".join(self.names) + "\n"
        _header_str += ",".join(self.types) + "\n"

        with open(self.filename,'w') as f:
            f.write(_header_str)

    def write_simulation_results(self,
            sim_id,
            results,
            filename=None):
        _sim_result = [str(sim_id)]
        _sim_result += [str(v) for k,v in results['parameters'].items()]
        _sim_result += [str(v) for k,v in results['qois'].items()]
        _sim_result += [str(v) for k,v in results['errors'].items()]

        _str_sim_results = ",".join(_sim_result)

        with open(filename,'a') as f:
            f.write(_str_sim_results+"\n")


    def read(self,filename=None):
        if filename is not None:
            self.filename = filename
        try:
            with open(self.filename,'r') as f:
                lines = f.readlines()
        except FileNotFoundError as e:
            print("cannot find file: {}".format(self.filename))
            print("current_working_dir: {}".format(os.getcwd()))
            raise

        self.names = [s.strip() for s in lines[0].strip().split(',')]
        self.types = [s.strip() for s in lines[1].strip().split(',')]

        table = []
        for line in lines[2:]:
            tokens = line.strip().split(',')
            values = []
            for i in range(len(self.names)):
                if i == 0:
                    values.append(tokens[i])
                else:
                    values.append(float(tokens[i]))
            table.append(values)
        
        self.df = pd.DataFrame(table)
        self.df.columns = self.names
        
        self.parameter_names = [
                n for i,n in enumerate(self.names) \
                    if self.types[i] == 'param']
        self.qoi_names = [
                n for i,n in enumerate(self.names) \
                        if self.types[i] == 'qoi']
        self.error_names = [
                n for i,n in enumerate(self.names) \
                        if self.types[i] == 'err']
        self.score_names = [
                n for i,n in enumerate(self.names) \
                        if self.types[i] == 'score']

        self.parameter_df = self.df[self.parameter_names]
        self.error_df = self.df[self.error_names]
        self.qoi_df = self.df[self.qoi_names]
        
    def write(self,filename):
        fn = filename

        s =  [",".join(self.names)]
        s += [",".join(self.types)]
        s += [",".join([str(v) for v in k]) for k in self.df.values.tolist()]

        with open (fn, 'w') as f:
            f.write("\n".join(s))

    def score_by_sum_if_less_than_median(
            self,
            error_df=None,
            err_type='abs'):

        _abs_error_df = self.error_df.copy(deep=True)
        _abs_error_df = _abs_error_df.loc[:,_abs_error_df.columns != 'sim_id'].abs()
        _sub_medians = _abs_error_df - _abs_error_df.median(axis=0)
        _scores = _sub_medians.copy(deep=True)
        _scores[_scores > 0] = 0 # no points if greater than median
        _scores[_scores < 0 ] = 1 # one point if less than median
        _metric = _scores.sum(axis=1)

        # calculate the metric
        self.names.append('sum_b_lt_median')
        self.score_names.append('sum_b_lt_median')
        self.types.append('score')
        self.df['sum_b_lt_median'] = np.copy(_metric)

    def score_by_d_metric(
            self,
            error_df=None,
            scaling_factors='DFT',
            err_type='abs'):

        _sf = 'd_metric'
        if type(scaling_factors) == str:
            _qoi_ref = scaling_factors
            if self.scaling_factors is None:
                self.scaling_factors = OrderedDict()
            self.scaling_factors[_sf] = OrderedDict()
            for qn in self.qoi_names:
                en = '{qoi_name}.err'.format(qoi_name=qn)
                self.scaling_factors[_sf][en] = 1/self.qoi_references[_qoi_ref][qn]
        elif isinstance(scaling_factors,dict):
            if self.scaling_factors is None:
                self.scaling_factors = OrderedDict()
            self.scaling_factors[_sf] = OrderedDict()
            for qn in self.qoi_names:
                self.scaling_factors[_sf][en] = scaling_factors[en]

        # normalize the errors
        _rescaled_error_df = None
        if err_type == 'abs':
            if error_df is None:
                _rescaled_error_df = self.error_df.copy(deep=True)
            else:
                _rescaled_error_df = error_df.copy(deep=True)
            for col in _rescaled_error_df:
                if col != 'sim_id':
                    sf = self.scaling_factors[_sf][col]
                    _rescaled_error_df[col] = sf*_rescaled_error_df[col].abs()
        self.rescaled_error_df = _rescaled_error_df.copy(deep=True)

        # calculate the metric
        _d_metric = np.sqrt(np.square(
                _rescaled_error_df.loc[:,_rescaled_error_df.columns != 'sim_id']
            ).sum(axis=1))
        self.names.append('d_metric')
        self.score_names.append('d_metric')
        self.types.append('score')
        self.df['d_metric'] = np.copy(_d_metric)

    def subselect_by_score(
            self,
            score_name,
            n=1000):
        """
            Args:
                scaling_factors(dict): the key is the error name, the value
                    is a scalar vaue from which the errors will be divided for
                    the purposes of scaling.
                n(int): the number of points to return
                result(str): should be either results, pareto or culled.
                    Default is culled.
        """

        # determine sub population
        if score_name == 'd_metric':
            self.sub_indices = self.df.nsmallest(n,'d_metric').index
        elif score_name == 'sum_b_lt_median':
            self.sub_indices = self.df.nlargest(n,'sum_b_lt_median').index
        else:
            err_msg = "{score_name} is not a valid score_name".format(
                    score_name=score_name)
            raise ValueError(err_msg)

        self.sub_df = self.df.loc[self.sub_indices]
        self.sub_error_df = self.error_df.loc[self.sub_indices]
        self.sub_parameter_df = self.parameter_df.loc[self.sub_indices]
        self.sub_qoi_df = self.qoi_df.loc[self.sub_indices]

    def write_subselect(self,filename=None):
        _filename = None
        if filename is None:
            _filename = "subselect.{}.out".format(
                    ".".join(self.score_names))
        elif type(filename) is str:
            _filename = filename
        else:
            err_msg = "the filename argument for this method must be a string"
            raise ValueError(err_msg)

        if type(self.sub_df) is None:
            err_msg="no subselection has been done on the data"
            raise ValueError(err_msg)
        if not isinstance(self.sub_df,pd.DataFrame) :
            print(type(self.sub_df))
            err_msg = "the sub_df attribute must be a pandas.DataFrame"
            raise ValueError(err_msg)

        # build the string
        str_out  = ','.join([n for n in self.names]) + "\n"
        str_out += ','.join([t for t in self.types]) + "\n"
        for row in self.sub_df.iterrows():
            _row = [a for i,a in enumerate(row[1])] #unpack tuple
            try:
                _row[0] = int(_row[0]) # row[0] is the sim_id
            except ValueError as e:
                raise
            #pass
            str_out += ','.join([str(s) for s in _row]) + "\n"

        with open(_filename,'w') as f:
            f.write(str_out)

        return _filename
