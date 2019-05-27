import pytest
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatDataAnalyzer

def testing_resources():
    testing_dict = OrderedDict()
    testing_dict['fn_config'] = 'resources/pyposmat.config.in'
    testing_dict['fn_data'] = 'resources/pyposmat.resutls.0.out'
    return testing_dict

def test____init____no_arguments():
    pda=PyposmatDataAnalyzer()
    assert type(pda) is PyposmatDataAnalyzer

def test____init____path_arguments():
    rscs=testing_resources()
    fn_config=rscs['fn_config']
    fn_data=rscs['fn_data']
    
    pda=PyposmatDataAnalyzer(
            fn_config=fn_config,
            fn_data=fn_data)

    assert type(pda) is PyposmatDataAnalyzer
    assert type(pda.configuration) is PyposmatConfigurationFile
    assert type(pda.datafile) is PyposmatDataFile
    assert type(pda.df) is pd.DataFrame
    assert type(pda.parameter_names) is list
    assert type(pda.error_names) is list
    assert type(pda.qoi_names) is list
    assert type(pda.names) is list

def test__read_configuration_file__validpath():
    rscs=testing_resources()
    fn_config=rscs['fn_config']

    pda=PyposmatDataAnalyzer()
    pda.read_configuration_file(filename=fn_config)

    assert type(pda.configuration) is PyposmatConfigurationfile

def test__read_data_file__validpath():
    rscs=testing_resources()
    fn_data=rscs['fn_data']

    pda=PyposmatDataAnalyzer()
    pda.read_data_file(filename=fn_data)

    assert type(pda.datafile) is PyposmatDataFile

def test__calculate_pareto_set(self,df=None,fn_out=None,sz=500):
    rscs=testing_resources()
    fn_data=rscs['fn_data']

    pda=PyposmatDataAnalyzer()
    pda.read_data_file(filename=fn_data)

    assert type(pda.datafile) is PyposmatDataFile
    
 


