import os
from collections import OrderedDict
import yaml
import pypospack.utils
from pypospack.io.filesystem import OrderedDictYAMLLoader
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatDataAnalyzer

class NewPyposmatDataAnalyzer(PyposmatDataAnalyzer):

    def str__analysis_results(self,analysis_results=None):

        assert analysis_results is None or isinstance(analysis_results,OrderedDict)

        if analysis_results is None:
            analysis_results = self.analysis_results

        assert isinstance(self.analysis_results,OrderedDict)

        s = ""
        for i in range(self.configuration.n_iterations):
            s += "iteration_results:{}/{}\n".format(i,self.configuration.n_iterations)
            s += "results_statistics:{}\n".format(analysis_results[i]['results_statistics'])
            s += "filter_info:{}\n".format(analysis_results[i]['filter_info'])
            s += "kde_statistics:{}\n".format(analysis_results[i]['kde_statistics'])

        return s

def get_testing_set():

    testing_set = OrderedDict()
    testing_set['data_directory'] = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'data','Si__sw__data','pareto_optimization_1'
    )
    testing_set['config_fn'] = os.path.join(
            testing_set['data_directory'],
           'pyposmat.config.in'
    )
    testing_set['analysis_fn'] = 'pyposmat.analysis.out'

    assert os.path.isdir(testing_set['data_directory'])
    assert os.path.isfile(testing_set['config_fn'])

    return testing_set

def dev__read_analysis_file():
    
    testing_set = get_testing_set()
    
    o = NewPyposmatDataAnalyzer()
    o.initialize_configuration(config_fn=testing_set['config_fn'])
    o.read_analysis_file(filename=testing_set['analysis_fn'])
    print(o.str__analysis_results())

def dev__write_kde_file():

    testing_set= get_testing_set()
    i_iteration = 1

    o = NewPyposmatDataAnalyzer()
    o.initialize_configuration(config_fn=testing_set['config_fn'])

def dev__analyze_kde_data():

    testing_set = get_testing_set()
    i_iteration = 1

    kde_fn = os.path.join(
            testing_set['data_directory'],
            'pyposmat.kde.{}.out'.format(i_iteration+1)
    )
    print('args:')
    print('\tkde_fn:{}'.format(kde_fn))
    print('\tos.path.isfile(kde_fn):{}'.format(os.path.isfile(kde_fn)))
        
    o = NewPyposmatDataAnalyzer()
    o.initialize_configuration(config_fn=testing_set['config_fn'])

    o.analyze_kde_data(
            i_iteration=i_iteration,
            filename=kde_fn
    )

    print("o.kde_statistics={}".format(o.kde_statistics))
    print(o.str__descriptive_statistics(o.kde_statistics))

def dev__str__analysis_results():
    testing_set = get_testing_set()

    o = NewPyposmatDataAnalyzer()
    o.initialize_configuration(config_fn=testing_set['config_fn'])
    o.read_analysis_file(filename=testing_set['analysis_fn'])

    print(o.str__analysis_results())
if __name__ == "__main__":
    #dev__read_analysis_file()
    dev__analyze_kde_data()
    dev__str__analysis_results()
