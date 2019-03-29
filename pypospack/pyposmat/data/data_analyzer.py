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


class PyposmatDataAnalyzer(object):
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
        self.kde_data = None
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
        self.descriptive_statistics = self.get_descriptive_statistics()

        self.filter_set_info = OrderedDict() 

        for filter_type, filter_info in self.configuration.qoi_constraints.items():
           
            if self.is_debug:
                print('filter_type:{}'.format(filter_type))
                print('filter_info:{}'.format(filter_info))

            if filter_type in ['qoi_constraints']:
                is_survive_idx, new_filter_info = self.filter_by_qoi_constraints()

                self.filter_set_info['filter_by_qoi_constraints'] = OrderedDict([
                    ('is_survive_idx', is_survive_idx),
                    ('filter_info', new_filter_info)
                ])

            elif filter_type in ['select_pareto_only','filter_by_pareto_membership']:
                is_survive_idx, new_filter_info = self.filter_by_pareto_membership()
                
                self.filter_set_info['filter_by_pareto_membership'] = OrderedDict([
                    ('is_survive_idx',is_survive_idx),
                    ('filter_info', new_filter_info)
                ])

            elif filter_type in ['filter_by__d_zerror','filter_by_d_zerror']:

                try:
                    n_potentials_min = filter_info['n_potentials_min']
                except KeyError as e:
                    n_potential_min = None

                try:
                    n_potentials_max = filter_info['n_potentials_max']
                except KeyError as e:
                    n_potentials_max = None

                is_survive_idx, new_filter_info = self.filter_by_cost_function(
                        loss_function_type='abs_error',
                        cost_function_type='weighted_sum',
                        weighting_scheme_type='scale_by_z_normalization',
                        pct_to_keep=filter_info['pct_to_keep'],
                        n_potentials_min=n_potentials_min,
                        n_potentials_max=n_potentials_max,
                )
             
                self.filter_set_info['filter_by_cost_function'] = OrderedDict([
                    ('is_survive_idx',is_survive_idx),
                    ('filter_info', new_filter_info)
                ])

            elif filter_type in ['filter_by_cost_function']:

                try:
                    n_potentials_min = filter_info['n_potentials_min']
                except KeyError as e:
                    n_potential_min = None

                try:
                    n_potentials_max = filter_info['n_potentials_max']
                except KeyError as e:
                    n_potentials_max = None
 
                is_survive_idx, new_filter_info = self.filter_by_cost_function(
                        loss_function_type=filter_info['loss_function_type'],
                        cost_function_type=filter_info['cost_function_type'],
                        weighting_scheme_type=filter_info['weighting_scheme_type'],
                        pct_to_keep=filter_info['pct_to_keep'],
                        n_potentials_min=n_potentials_min,
                        n_potentials_max=n_potentials_max
                )

                self.filter_set_info['filter_by_cost_function'] = OrderedDict([
                    ('is_survive_idx',is_survive_idx),
                    ('filter_info', new_filter_info)
                ])

            else:
                m =  "unknown filtering type to constraint candidate potentials."
                m += "filter_type={}".format(filter_type)
                raise PyposmatUnknownQoiFilterType(m)

        all_is_survive_idx = [v['is_survive_idx'] for v in self.filter_set_info.values()]
        self.filter_set_info['is_survive_idx'] = set.intersection(*all_is_survive_idx)

        n_potentials_start = self.results_df.shape[0]
        n_potentials_final = len(self.filter_set_info['is_survive_idx'])
        pct_to_keep = n_potentials_final/n_potentials_start
        self.filter_set_info['n_potentials_start'] = n_potentials_start
        self.filter_set_info['n_potentials_final'] = n_potentials_final
        self.filter_set_info['pct_to_keep'] = pct_to_keep

        return self.filter_set_info['is_survive_idx'], self.filter_set_info
   
    def write_kde_file(self,filename):
        names = ['sim_id'] + self.parameter_names + self.qoi_names + self.error_names
       
        self.kde_data = PyposmatDataFile()
        self.kde_data.read(filename=self.results_data_fn)
        self.kde_data.df = self.kde_data.df.iloc[list(self.filter_set_info['is_survive_idx'])]
        self.kde_data.write(filename=filename)


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

    def get_weights_by_qoi_target_normalization(self,df=None):
        """ determine weights using qoi target normalization

        Args:
            df(pandas.DataFrame,None):  The dataframe constraining the data.  The default
                value is Noen, which sets the arguement value to the property results_df.
        """

        assert df is None or isinstance(pd.DataFrame)

        if df is None:
            df = self.results_df

        weights = OrderedDict([(k,abs(1/v)) for k,v in self.configuration.qoi_targets.items()])
        sum_weights = sum([v for v in weights.values()])
        for qoi_name, qoi_weight in weights.items():
            weights[qoi_name] = qoi_weight/sum_weights

        return weights


    def filter_by_cost_function(self,
                                df=None,
                                loss_function_type='abs_error',
                                cost_function_type='weighted_sum',
                                weighting_scheme_type='scale_by_qoi_target',
                                pct_to_keep=None,
                                n_potentials_min=50,
                                n_potentials_max=10000):
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
        assert weighting_scheme_type in self.weighting_scheme_types
        assert isinstance(pct_to_keep,float) and pct_to_keep <= 1

        if df is None:
            df = self.results_df

        weighting_to_method_map = OrderedDict([
            ('scale_by_qoi_target',self.get_weights_by_qoi_target_normalization),
            ('scale_by_z_normalization',self.get_weights_by_z_error_normalization)
        ])

        loss_function_to_method_map = OrderedDict([
            ('abs_error',self.create_absolute_errors),
            ('squared_error',self.create_squared_errors)
        ])

        cost_function_to_method_map = OrderedDict([
            ('weighted_sum',self.get_weighted_cost_function)
        ])

        try:
            weights = weighting_to_method_map[weighting_scheme_type]()
        except KeyError as e:
            m = "{}: unknown weighting scheme".format(weighting_scheme_type)
            raise

        assert all([k in self.qoi_names for k in weights])
        assert all([isinstance(v,float) for v in weights.values()])

        try:
            loss_error_names, loss_function_df = loss_function_to_method_map[loss_function_type]()
        except KeyError as e:
            m = "{}: unknown loss function type".format(loss_error_names)
            raise

        try:
            df['cost'] = cost_function_to_method_map[cost_function_type](
                    df=loss_function_df,
                    weights=weights,
                    loss_error_names=loss_error_names
            )

        except KeyError as e:
            m = '{}: unknown cost function type'.format(loss_error_names)
            raise

        n_potentials_start, n_cols = df.shape
        n_potentials_final = int(n_potentials_start*pct_to_keep)
        n_potentials_final = min(n_potentials_start,max(n_potentials_min,n_potentials_final))
        n_potentials_final = min(n_potentials_max,n_potentials_final)
        is_survive_idx = df['cost'].nsmallest(n_potentials_final).index
        is_survive_idx = is_survive_idx.values.tolist()
        is_survive_idx.sort()

        filter_info = OrderedDict()
        filter_info['filter_type'] = 'filter_by_cost_function'
        filter_info['loss_function_type'] = loss_function_type
        filter_info['cost_function_type'] = cost_function_type
        filter_info['weights'] = weights
        filter_info['pct_to_keep'] = pct_to_keep
        filter_info['n_potentials_min'] = n_potentials_min
        filter_info['n_potentials_max'] = n_potentials_max
        filter_info['n_potentials_start'] = n_potentials_start
        filter_info['n_potentials_final'] = n_potentials_final

        return set(is_survive_idx),filter_info

    def create_absolute_errors(self,df=None): 
        
        if df is None:
            df = self.results_df

        loss_error_names = ['{}.abserr'.format(v) for v in self.qoi_names]
        error_names = ['{}.err'.format(v) for v in self.qoi_names]

        df[loss_error_names] = df[error_names].abs()

        return loss_error_names, df[loss_error_names]

    def create_squared_errors(self,df=None):

        if df is None:
            df = self.results_df

        loss_error_names = ['{}.sqerr'.format(v) for v in self.qoi_names]
        error_names = ['{}.err'.format(v) for v in self.qoi_names]

        df[loss_error_names] = df[error_names].square()

        return loss_error_names, df[loss_error_names]

    def get_weighted_cost_function(self,weights,loss_error_names,df=None):

        if df is None:
            df = self.results_df

        weights_df = pd.DataFrame(
                pd.Series([v for v in weights.values()],
                          index=['{}.abserr'.format(k) for k in weights],
                          name=0
                )
        )

        cost_df = df[loss_error_names].dot(weights_df)
       
        return cost_df

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

        if len(is_survive_idx) > 0:
            is_survive_idx = set.intersection(*[set(v) for v in is_survive_idx])
            n_potentials_final = len(is_survive_idx)
        else:
            n_potentials_final = n_potentials_start

        qoi_constraint_info['n_potentials_final'] = n_potentials_final
        qoi_constraint_info['pct_potentials_final'] = n_potentials_final/n_potentials_start

        return is_survive_idx,qoi_constraint_info

    def str__analysis_results(self):
        
        filter_to_str_method_map = OrderedDict([
            ('filter_by_pareto_membership',self.str__filter_by_pareto_membership),
            ('filter_by_cost_function',self.str__filter_by_cost_function),
            ('filter_by_qoi_constraints',self.str__filter_by_qoi_constraints)
        ])

        s = ""
        # formatted output for descriptive statistics
        s += '{}\n'.format(self.str__descriptive_statistics(
            descriptive_statistics=self.descriptive_statistics
        ))

        # formatted output for filtering information
        s += 80*'-' + "\n"
        s += '{:^80}\n'.format('qoi filtering summary')

        keys_to_skip = [
                'n_potentials_start','n_potentials_final','is_survive_idx','pct_to_keep'
        ]
        for k,v in self.filter_set_info.items():
            if k in keys_to_skip:
                pass
            else:
                try:
                    s += '{}\n'.format(filter_to_str_method_map[k](v['filter_info']))
                except KeyError as e:
                    s += '{} {}\n'.format(k,v)
                    raise

        # formatted output for final results
        s += 80*'-' + "\n"
        s += 'final analysis\n'
        s += len('final_analysis')*'-' + "\n"
        s += "n_potentials_start:{}\n".format(self.filter_set_info['n_potentials_start'])
        s += "n_potentials_final:{}\n".format(self.filter_set_info['n_potentials_final'])
        s += "pct_to_keep:{}\n".format(self.filter_set_info['pct_to_keep'])

        return(s)

    def str__filter_by_pareto_membership(self,filter_info):
        s = 80*'-' + "\n"
        s += "{:<80}\n".format('filter_by_pareto_membership')
        s += len('filter_by_pareto_membership')*'-' + "\n"
        for k,v in filter_info.items():
            s += "{}:{}\n".format(k,v)

        return s

    def str__cost_function_weights(self,weights):
        
        assert isinstance(weights,OrderedDict)

        s = "{:^20} {:^10}\n".format('qoi_name','qoi_weight')
        s += "{} {}\n".format(20*'-',10*'-')
        for qoi_name, qoi_weight in weights.items():
            s+= "{:^20} {:^10.4e}\n".format(qoi_name,qoi_weight)
        return s

    def str__filter_by_cost_function(self,filter_info):
        is_debug = False

        if is_debug:
            for k,v in filter_info.items():
                print('{}:{}'.format(k,v))

        s = 80*'-' + "\n"
        s += "{:<80}\n".format('filter_by_cost_function')
        s += len('filter_by_cost_function')*'-' + "\n"

        s += 'loss_function_type:{}\n'.format(filter_info['loss_function_type'])
        s += 'cost_function_type:{}\n'.format(filter_info['cost_function_type'])
        
        s += self.str__cost_function_weights(weights=filter_info['weights'])

        s += "pct_to_keep:{}\n".format(filter_info['pct_to_keep'])
        s += "n_potentials_min:{}\n".format(filter_info['n_potentials_min'])
        s += "n_potentials_max:{}\n".format(filter_info['n_potentials_max'])
        s += "n_potentials_start:{}\n".format(filter_info['n_potentials_start'])
        s += "n_potentials_final:{}\n".format(filter_info['n_potentials_final'])

        return s

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

    def str__filter_by_qoi_constraints(self,qoi_constraint_info):
        s = 80*'-' + "\n"
        s += "{:<80}\n".format('filter_by_qoi_constraints')
        s += len('filter_by_qoi_constraints')*'-' + "\n"

        s += '{:20} {:10} {:10}\n'.format('short_description','n_survived','pct_survived')
        for k,v in qoi_constraint_info['constraints'].items():
            short_description = v['short_description']
            n_survived = v['n_survived']
            pct_survived = v['pct_survived']
            s += '{:20} {:10} {:10.4%}\n'.format(short_description,n_survived,pct_survived)
        s += 'n_potentials_start={}\n'.format(qoi_constraint_info['n_potentials_start'])
        s += "n_potentials_final={}\n".format(qoi_constraint_info['n_potentials_final'])
        s += "pct_potentials_final={}\n".format(qoi_constraint_info['pct_potentials_final'])
        return s
         
