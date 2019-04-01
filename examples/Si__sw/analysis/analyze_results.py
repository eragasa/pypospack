import os
from collections import OrderedDict
import yaml
import pypospack.utils
from pypospack.io.filesystem import OrderedDictYAMLLoader
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatDataAnalyzer

data_directory = os.path.join(
        pypospack.utils.get_pypospack_root_directory(),
        'data','Si__sw__data','pareto_optimization_1'
)
config_fn = os.path.join(data_directory,'pyposmat.config.in')
analysis_fn = 'pyposmat.analysis.out'

assert os.path.isdir(data_directory)
assert os.path.isfile(config_fn)

o = PyposmatDataAnalyzer()
o.initialize_configuration(config_fn=config_fn)

n_iterations = o.configuration.n_iterations
print("n_iterations:{}".format(n_iterations))

# read analysis results to file system
if os.path.isfile(analysis_fn):
    with open(analysis_fn,'r') as f:
        analysis_results = yaml.load(f, OrderedDictYAMLLoader)
else:
    analysis_results = OrderedDict()

for i in range(n_iterations):
    if i not in analysis_results:
        print('i_iteration:{}'.format(i))
        analysis_results[i] = OrderedDict()
        analysis_results[i]['results_statistics'] = None
        analysis_results[i]['kde_statistics'] = None
        analysis_results[i]['filter_info'] = None

        filter_info = o.analyze_results_data(
                i_iteration=i,
                filename=os.path.join(
                    data_directory,'pyposmat.results.{}.out'.format(i)
                )
        )
        o.analyze_kde_data(
            i_iteration=i,
            filename=os.path.join(
                data_directory,'pyposmat.kde.{}.out'.format(i+1)
            )
        )

        analysis_results[i]['results_statistics'] = o.results_statistics
        analysis_results[i]['kde_statistics'] = o.kde_statistics
        analysis_results[i]['filter_info'] = filter_info

        print(analysis_results[i]['results_statistics'])
        print(analysis_results[i]['filter_info'])
        print(analysis_results[i]['kde_statistics'])
    else:
        pass

# write analysis results to filesystem
with open(analysis_fn,'w') as f:
    yaml.dump(analysis_results, f, default_flow_style=False)

for i in range(n_iterations):
    print("{}".format(i))
    for k,v in analysis_results[i].items():
        print("\t{}".format(k))
if __name__ == "__main__":
    pass
