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

from pypospack.exceptions import PyposmatUnknownQoiFilterType
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

class PyposmatAnalysisResults(object):

    def __init__(self,config_fn=None,o_config=None):

        assert config_fn is None or isintance(config_fn,str)
        assert o_config is None or isintance(o_config,PyposmatConfigurationFile)

        self.is_debug = False
        self.config_fn = None
        self.configuration = None
        self.iteration_results = OrderedDict()

        self.initialize_configuration(config_fn=config_fn,o_config=o_config)

    def read(self,filename):
        pass

    def write(self,filename):
        pass

    def initialize_configuration(self,config_fn=None,o_config=None):
        """ read the configuration file

        Must provide either the path of the configuration file through the config_fn
        argument, or provide an instance of a PyposmatConfigurationFile object.  If no
        arguments are provided, it will read the config_fn attribute as the path of the
        configuration file.

        Args:
            config_fn(str,None): the path of the configuration file to read.  By default,
               this attribute is None.
            o_config(str,None): an instance of a PyposmatConfigurationFile object.  By default,
               this attribute is None.

        """

        assert config_fn is None or isinstance(config_fn,str)
        assert o_config is None or isinstance(o_config,PyposmatConfigurationFile)

        if config_fn is not None and o_config is not None:
            m = (
                "must either provide the path to config_fn or a PyposmatConfigurationFile "
                "instance to to o_config"
            )
            raise TypeError(m)
        # default behavior
        elif (config_fn is None) and (o_config is None):
            self.configuration = PyposmatConfigurationFile()
            self.configuration.read(filename=self.config_fn)
        # a path is provided
        elif isinstance(config_fn,str):
            self.config_fn = config_fn
            self.configuration = PyposmatConfigurationFile()
            self.configuration.read(filename=self.config_fn)
        # an object is provided
        elif isinstance(o_config,PyposmatConfigurationFile):
            self.config_fn = None
            self.configuration = o_config
        else:
            m = (
                "must either provide the path to config_fn or a PyposmatConfigurationFile "
                "instance to to o_config"
            )
            raise TypeError(m)

class NewPyposmatDataAnalyzer(object):
    """ class to analyze the results of the simulations 
    
    Args:
        config_fn(str): path to the configuration file
        results_data_fn(str): path to the data file in which to do 
            analysis

    Attributes:
        config_fn(str): path to the configuration file
        results_fn(str): path to the results data file
        analysis_fn(str): path to the analysis data file
        configuration(PyposmatConfigurationFile): instance of the configuration
            file
        results_data(PyposmatDataFile): instance of the configuration file
        analysis(PyposmatAnalysisFile): instanceo of the analysis file which is used to marshal 
            information from previous iterations, and unmarshall information to file for later
            analysis
        mpi_rank(int): the MPI rank of the process running this code
        mpi_size(int): the MPI rank of the process running this code
    """
    
    loss_function_types = ['abs_error','squared_error']
    cost_function_types = ['weighted_sum','weighted_sum_product']
    weighting_scheme_types = ['scale_by_qoi_target','scale_by_z_normalization']
    
    def __init__(self,
                 config_fn=None,
                 results_data_fn=None,
                 analysis_fn = None,
                 o_config=None,
                 o_results_data=None,
                 mpi_rank=0,
                 mpi_size=1):

        # check arguments are the correct type
        assert config_fn is None or isinstance(config_fn,str)
        assert results_data_fn is None or isinstance(results_data_fn,str)
        assert analysis_fn is None or isinstance(analysis_fn,str)
        assert o_config is None or isinstance(o_config,PyposmatConfigurationFile)
        assert o_results_data is None or isinstance(o_results_data,PyposmatDataFile)

        # public attributes are set to NoneType objects for debugging purposes
        self.configuration = None
        self.results_data = None
        self.config_fn = None
        self.results_data_fn = None
        self.analysis_fn = None
        self.analysis = None
     
        # set debug level
        self.is_debug = False
        
        # set mpi info
        self.mpi_rank = mpi_rank
        self.mpi_size = mpi_size

        # configuration file
        if isinstance(config_fn,str) or isinstance(o_config,PyposmatConfigurationFile):
            self.config_fn = config_fn
            self.initialize_configuration(config_fn=config_fn,
                                          o_config=o_config
            )
        
        # results data file
        if isinstance(results_data_fn,str) or isinstance(o_results_data,PyposmatDataFile):
            self.results_data_fn = results_data_fn
            self.initialize_results_data(results_data_fn=results_data_fn,
                                         o_results_data=o_results_data
            )

    @property
    def parameter_names(self):
        if isinstance(self.configuration,PyposmatConfigurationFile):
            return self.configuration.parameter_names
        else:
            return None

    @property
    def error_names(self):
        if isinstance(self.configuration,PyposmatConfigurationFile):
            return self.configuration.error_names
        else:
            return None

    @property
    def qoi_names(self):
        if isinstance(self.configuration,PyposmatConfigurationFile):
            return self.configuration.qoi_names
        else:
            return None

    @property
    def results_df(self):
        if isinstance(self.results_data,PyposmatDataFile):
            return self.results_data.df
        else:
            return None

    @results_df.setter
    def results_df(self,df):
        assert isinstance(df,pd.DataFrame)
        assert isintance(self.results_data,PyposmatDataFile)
        self.results_data.df = df

    @property
    def n_potentials_start(self):
        if isinstance(self.results_data,PyposmatDataFile):
            return self.results_df.shape[0]
        else:
            return None

    def initialize_configuration(self,config_fn=None,o_config=None):
        """ read the configuration file

        Must provide either the path of the configuration file through the config_fn
        argument, or provide an instance of a PyposmatConfigurationFile object.  If no
        arguments are provided, it will read the config_fn attribute as the path of the
        configuration file.

        Args:
            config_fn(str,None): the path of the configuration file to read.  By default,
               this attribute is None.
            o_config(str,None): an instance of a PyposmatConfigurationFile object.  By default,
               this attribute is None.

        """

        assert config_fn is None or isinstance(config_fn,str)
        assert o_config is None or isinstance(o_config,PyposmatConfigurationFile)

        if config_fn is not None and o_config is not None:
            m = (
                "must either provide the path to config_fn or a PyposmatConfigurationFile "
                "instance to to o_config"
            )
            raise TypeError(m)
        # default behavior
        elif (config_fn is None) and (o_config is None):
            self.configuration = PyposmatConfigurationFile()
            self.configuration.read(filename=self.config_fn)
        # a path is provided
        elif isinstance(config_fn,str):
            self.config_fn = config_fn
            self.configuration = PyposmatConfigurationFile()
            self.configuration.read(filename=self.config_fn)
        # an object is provided
        elif isinstance(o_config,PyposmatConfigurationFile):
            self.config_fn = None
            self.configuration = o_config
        else:
            m = (
                "must either provide the path to config_fn or a PyposmatConfigurationFile "
                "instance to to o_config"
            )
            raise TypeError(m)

    def initialize_results_data(self,results_data_fn=None,o_results_data=None):
        """ read results data file

        This initializes the results_data attribute as PyposmatDataFile object.  There are 
        three mechanisms to initialize this object.  If default None arguments, are passed
        into this method, it will use the results_fn attribute.  Otherwise, it will use
        the not None arguement.  If two arguments are provided, a TypeError will be raised.

        Args:
            results_data_fn(str,None): the path to the results file.  By default it is None.
            o_results_data(PyposmatDataFile,None): the path to the results data file.  By 
                default it is None
        """

        if results_data_fn is not None and o_results_data is not None:
            m = (
                "must either provide the path to results_data_fn or a PyposmatDataFile"
                "instance to to o_results_data"
            )
            raise TypeError(m)
        # default behavior
        elif results_data_fn is None and o_results_data is None:
            self.results_data = PyposmatDataFile()
            self.results_data.read(filename=self.results_data_fn)
        # a path is provided
        elif isinstance(results_data_fn,str):
            self.results_data_fn = results_data_fn
            self.results_data = PyposmatDataFile()
            self.results_data.read(filename=self.results_data_fn)
        # an object is provided
        elif isinstance(o_results_data,PyposmatDataFile):
            self.results_data_fn = None
            self.results_data = o_results_data
        else:
            m = (
                "must either provide the path to results_data_fn or a PyposmatDataFile"
                "instance to to o_results_data"
            )
            raise TypeError(m)

    @property
    def qoi_constraints(self):
        if isinstance(self.configuration,PyposmatConfigurationFile):
            if 'qoi_constraints' in self.configuration.qoi_constraints:
                return self.configuration.qoi_constraints['qoi_constraints']
            else:
                return None
        else:
            return None
   
    def read_analysis_file(self,filename=None):
        """ read the pypospack analysis file

        Args:
            filename(str,None): By default this is set to None, which uses the attribute
               analysis_fn

        """
        if isinstance(filename,str):
            self.analysis_fn = filename
            self.analysis = PyposmatAnalysisFile(o_config=self.configuration)
            self.analysis.read(filename=filename)

    def get_descriptive_statistics(self,df=None):

        assert df is None or isinstance(df,pd.DataFrame)

        if df is None:
            df = self.results_df
        descriptive_statistics =OrderedDict()
        descriptive_statistics['qois'] = OrderedDict()
        
        for qn in self.qoi_names:
            descriptive_statistics['qois'][qn] = OrderedDict()
            descriptive_statistics['qois'][qn]['target'] = self.configuration.qoi_targets[qn]
            descriptive_statistics['qois'][qn]['mean'] = df[qn].mean()
            descriptive_statistics['qois'][qn]['std'] = df[qn].std()

        return descriptive_statistics

    def str__descriptive_statistics(self,descriptive_statistics):
        """ convert desciprtive statistics to a display to string

        Args:
            descriptive_statistics(OrderedDict): a nested ordered dictionary containing
                the descriptive statistics of the last simulation
        Rerturn:
            str
        """
        
        len_qoi_names = '{:20}'
        len_qoi_target = '{:+10.4e}'
        len_qoi_mean = '{:+10.4e}'
        len_qoi_std = '{:+10.4e}'
        QOI_HEADER_FORMAT = '{:^20} {:^11} {:^11} {:^11}\n'

        QOI_ROW_FORMAT = " ".join([len_qoi_names,len_qoi_target,len_qoi_mean,len_qoi_std])+"\n"

        s = 80*"-" + "\n"
        s += "{:^80}\n".format('Descriptive Statistics')
        s += 80*"-" + "\n"
        s += "\n"
        s += QOI_HEADER_FORMAT.format('qoi_names','qoi_target','mean','std')
        for k,v in descriptive_statistics['qois'].items():
            s += QOI_ROW_FORMAT.format(k,v['target'],v['mean'],v['std'])
        return s

    def analyze_results(self,i_iteration,filename=None,o_data=None):

        assert isinstance(i_iteration,int)
        assert filename is None or isinstance(filename,str)
        assert o_data is None or isinstance(o_data, PyposmatDataFile)

        if isinstance(filename,str) or isinstance(o_data,PyposmatDataFile):
            self.initialize_data_file(results_data_fn=filename,o_results_data=o_data)

        assert isinstance(self.results_data,PyposmatDataFile)
        assert i_iteration < self.configuration.n_iterations

        if False:
            # don't have testing files to analyze of multiple iterations
            sim_ids = self.results_df['sim_id'].values.tolist()
            print(sim_ids)

        descriptive_statistics = self.get_descriptive_statistics()

        filter_set_info = OrderedDict() 
        for filter_type,filter_info in self.configuration.qoi_constraints.items():
           
            if self.is_debug:
                print('filter_type:{}'.format(filter_type))
                print('filter_info:{}'.format(filter_info))
            if filter_type in ['qoi_constraints']:
                is_survive_idx, filter_info = self.filter_by_qoi_constraints()
                filter_set_info['filter_by_qoi_constraints'] = OrderedDict()
                filter_set_info['filter_by_qoi_constraints']['is_survive_idx'] = is_survive_idx
                filter_set_info['filter_by_qoi_constraints']['filter_info'] = qoi_constraint_info

            elif filter_type in ['select_pareto_only','filter_by_pareto_membership']:
                is_survive_idx, filter_info = self.filter_by_pareto_membership()
                filter_set_info['filter_by_pareto_membership'] = OrderedDict()
                filter_set_info['filter_by_pareto_membership']['is_survive_idx'] = is_survive_idx
                filter_set_info['filter_by_pareto_membership']['filter_info'] = pareto_set_info

            elif filter_type in ['filter_by__d_zerror','filter_by_d_zerror']:
                pct_to_keep = filter_info['percentile']
                weights = self.get_weights_by_z_error_normalization()

                is_survive_idx,cost_function_info = self.filter_by_cost_function(
                        loss_function_type='abs_error',
                        cost_function_type='weighted_sum',
                        weights=weights,
                        pct_to_keep=pct_to_keep
                )
                
                temp_filter_info = OrderedDict()
                temp_filter_info['loss_function_type'] = 'abs_error'
                temp_filter_info['cost_function_type'] = 'weighted_sum'
                temp_filter_info['weighting_scheme'] = 'scale_by_z_normalization'
                temp_filter_info['weights'] = weights
                temp_filter_info['pct_to_keep'] = pct_to_keep
                filter_set_info['filter_by_cost_function'] = OrderedDict()
                filter_set_info['filter_by_cost_function']['is_survive_idx'] = is_survive_idx
                filter_set_info['filter_by_cost_function']['filter_info'] = temp_filter_info
                filter_set_info['filter_by_cost_function']['filter_info']['loss_function_type']

            elif filter_type in ['filter_by_cost_function']:

                # get weights based on the weighting scheme
                weighting_scheme_type = filter_info['weighting_scheme_type']
                assert weighting_scheme_type in self.weighting_scheme_types

                weighting_to_method_map = OrderedDict([
                    ('scale_by_qoi_target',self.get_weights_by_qoi_target_normalization),
                    ('scale_by_z_normalization',self.get_weights_by_z_error_normalization)
                ])

                try:
                    weights = weighting_to_method_map[weighting_scheme_type]()
                except KeyError as e:
                    m = "{}: unknown weighting scheme".format(weighting_scheme_type)
                    raise
               
                print(weights)
                exit()
                
                if filter_info['cost_function_type'] == 'weighted_sum':
                    pass

                if filter_info['loss_function_type'] == 'absolute_error':
                    pass
                elif filter_info['loss_function_type'] == 'squared_error':
                    pass

            else:
                m =  "unknown filtering type to constraint candidate potentials."
                m += "filter_type={}".format(filter_type)
                raise PyposmatUnknownQoiFilterType(m)

        all_is_survive_idx = [v['is_survive_idx'] for v in filter_set_info.values()]
        filter_set_info['is_survive_idx'] = set.intersection(*all_is_survive_idx)
        filter_set_info['n_potentials_start'] = n_potentials_start
        filter_set_info['n_potentials_final'] = len(filter_set_info['is_survive_idx'])

        return filter_set_info['is_survive_idx'], filter_set_info
    
    def str__filter_by_cost_function_filter_info(self,filter_info):
        s = ""
        for k,v in filter_info.items():
            if k not in ['n_potentials_start','n_potentials_end','weights']:
                s += "{}:{}\n".format(k,v)
            if k == 'weights':
                s += "{:^20} {:^20}\n".format('qoi_name','qoi_weight')
                s += "{} {}\n".format(20*'-',20*'-')
                for qoi_name, qoi_weight in v.items():
                    s += "{:^20} {:20.4}\n".format(qoi_name,qoi_weight)
            else:
                s += "{}:{}\n".format(k,v)
        return s

    def get_weights_by_z_error_normalization(self,df=None):
        """ determine weights using z normalization

        Args:
            df(pandas.DataFrame,None): The dataframe constaining the data.  The default
                value is None, which sets the argument value to the property results_df.
        """

        assert df is None or isinstance(pd.DataFrame)

        if df is None:
            df = self.results_df

        weights = OrderedDict([(k,1/df[k].std()) for k in self.qoi_names])
        
        return weights


    def filter_by_cost_function(self,
                                df=None,
                                loss_function_type='abs_error',
                                cost_function_type='weighted_sum',
                                weights=None,
                                pct_to_keep=None):
        """
        Args:
            df(pandas,DataFrame): a pandas dataframe object contain the results of the this
                iterations of sampling
            loss_function_type(str): the type of loss function to use
            cost_function_type(str): the type of cost function to use
        """
 
        assert df is None or isinstance(df,pd.DataFrame)
        assert loss_function_type in self.loss_function_types
        assert cost_function_type in self.cost_function_types
        assert isinstance(weights,OrderedDict)
        assert all([k in self.qoi_names for k in weights])
        assert all([isinstance(v,float) for v in weights.values()])
        assert isinstance(pct_to_keep,float) and pct_to_keep <= 1

        if df is None:
            df = self.results_df


        if loss_function_type == 'abs_error':
            loss_error_names = ['{}.abserr'.format(v) for v in self.qoi_names]
            for v in self.qoi_names:
                df['{}.abserr'.format(v)] = df['{}.err'.format(v)].abs()
            weights_df = pd.DataFrame(
                    pd.Series([v for v in weights.values()],
                              index=['{}.abserr'.format(k) for k in weights],
                              name=0
                    )
            )
        elif loss_function_type == "squared_name":
            loss_error_names = ['{}.sqerr'.format(v) for v in self.qoi_names]
            for v in self.qoi_names:
                df['{}.sqerr'.format(v)] = df['{}.err'.format(v)].square()
        else:
             s = "{} is an unsupported loss_function_type".format(loss_function_type)
             raise ValueError(m)

        if cost_function_type == 'weighted_sum':
            df['cost'] = df[loss_error_names].dot(weights_df)

        n_potentials_begin, n_cols = df.shape
        n_potentials_final = int(n_potentials_begin*pct_to_keep)
        is_survive_idx = df['cost'].nsmallest(n_potentials_final).index
        is_survive_idx = is_survive_idx.values.tolist()
        is_survive_idx.sort()

        filter_info = OrderedDict()
        filter_info['filter_type'] = 'filter_by_cost_function'
        filter_info['loss_function_type'] = loss_function_type
        filter_info['cost_function_type'] = cost_function_type
        filter_info['weights'] = weights
        filter_info['n_potentials_begin'] = n_potentials_begin
        filter_info['n_potentials_final'] = n_potentials_final
        return set(is_survive_idx),filter_info
        
    def filter_by_pareto_membership(self,df=None,v=None):
        """ filter by pareto membership

        Calculation of pareto membership by removal of dominated points
        Args:
            df(pandas.DataFrame): the dataframe of results.  By default this is set to None,
                which uses the results_df property from a PyposmatDataFile which was
                read into this instance.

        Returns:
            (list): a list of integers whic hare the indices of the of the surviving
                candidate parameterization
            (OrderedDict): a dictionary containing the information for Pareto membership
        """

        assert df is None or isintance(df,pd.DataFrame)

        if df is None:
            df = self.results_df
    
        for qn in self.qoi_names:
            df['{}.abserr'.format(qn)] = df['{}.err'.format(qn)].abs()
      
        abs_error_names = ['{}.abserr'.format(v) for v in self.qoi_names]
        is_survive_idx = pareto(df[abs_error_names].values.tolist())
        
        n_potentials_start, n_cols = df.shape
        n_potentials_final = len(is_survive_idx)
        pct_potentials_final = n_potentials_final/n_potentials_start
        pareto_set_info = OrderedDict()
        pareto_set_info['n_potentials_start'] = n_potentials_start
        pareto_set_info['n_potentials_final'] = n_potentials_final
        pareto_set_info['pct_potentials_final'] = pct_potentials_final

        return set(is_survive_idx), pareto_set_info

    def filter_by_qoi_constraints(self,df=None):
        """filter by qoi constraints

        Args:
            df(pd.DataFrame,None): 
        """

        if df is None:
            df = self.results_data.df

        for qn in self.qoi_names:
            df['{}.abserr'.format(qn)] = df['{}.err'.format(qn)].abs()

        n_potentials_start, n_cols = df.shape
        qoi_constraint_info= OrderedDict()
        qoi_constraint_info['n_potentials_start'] = n_potentials_start
        qoi_constraint_info['n_potentials_final'] = None
        qoi_constraint_info['pct_potentials_final'] = None
        qoi_constraint_info['constraints'] = OrderedDict()
        is_survive_idx = []
        for k,v in self.qoi_constraints.items():
            if v[0] == '<':
                is_survive_idx.append(df.index[df[k] < v[1]])
            elif v[0] == '>':
                is_survive_idx.append(df.index[df[k] > v[1]])
            elif v[0] == '=':
                is_survive_idx.append(df.index[df[k] > v[1]])
            else:
                raise ValueError('unknown operator for {} {} {}'.format(k,v[0],v[1]))
            n_survived = len(is_survive_idx[len(is_survive_idx)-1])

            constraint_info = OrderedDict()
            constraint_info['short_description'] = '{} {} {}'.format(k,v[0],v[1])
            constraint_info['n_survived'] = n_survived
            constraint_info['pct_survived'] = n_survived/n_potentials_start
            qoi_constraint_info['constraints'][k] = constraint_info

        is_survive_idx = set.intersection(*[set(v) for v in is_survive_idx])
        n_potentials_final = len(is_survive_idx)
        qoi_constraint_info['n_potentials_final'] = n_potentials_final
        qoi_constraint_info['pct_potentials_final'] = n_potentials_final/n_potentials_start

        return is_survive_idx,qoi_constraint_info

    def str__qoi_constraint_info(self,qoi_constraint_info):
        n_potentials_start = qoi_constraint_info['n_potentials_start']
        s = 'n_potentials_start={}\n'.format(n_potentials_start)
        for k,v in qoi_constraint_info['constraints'].items():
            short_description = v['short_description']
            n_survived = v['n_survived']
            pct_survived = n_survived/n_potentials_start
            s += '{} {} {}\n'.format(short_description,n_survived,pct_survived)
        s += "n_potentials_final={}\n".format(qoi_constraint_info['n_potentials_final'])
        s += "pct_potentials_final={}\n".format(qoi_constraint_info['pct_potentials_final'])
        return s
         
        
        

class PyposmatDataAnalyzer(object):
    """ class to analyze the results of the simulations.

    """
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
        return OrderedDict([(qn,qv['target'] )for qn,qv in self.configuration.qois.items()])

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
