import os
import copy
from collections import OrderedDict
import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.visualization.parallel_plot_new import PyposmatParallelCoordinatesPlot

parallel_plot_config= OrderedDict()
parallel_plot_config['args'] = OrderedDict()
parallel_plot_config['p_3.5_q_0.5'] = OrderedDict()
parallel_plot_config['p_3.5_q_0.5']['label'] = 'p=3.5 q=0.5'
parallel_plot_config['p_3.5_q_0.5']['color'] = 'blue'
parallel_plot_config['p_3.5_q_0.5']['data_directory'] = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),'data','Si__sw__data',
            'pareto_optimization_p_3.5_q_0.5')
parallel_plot_config['p_3.5_q_0.5']['config_fn'] = os.path.join(
            parallel_plot_config['p_3.5_q_0.5']['data_directory'],
            'pyposmat.config.in')
parallel_plot_config['p_3.5_q_0.5']['data_fn'] = os.path.join(
            parallel_plot_config['p_3.5_q_0.5']['data_directory'],
            'pyposmat.kde.19.out')

if __name__ == "__main__":
    # initialization
    o_plot = PyposmatParallelCoordinatesPlot()

    # add data to plot
    for k,v in parallel_plot_config.items():
        print(k,v)
        if k == 'args':
            pass
        else:
            o_config = PyposmatConfigurationFile()
            o_config.read(filename=v['config_fn'])
            
            o_data = PyposmatDataFile()
            o_data.read(filename=v['data_fn'])
            o_data.create_normalized_errors(qoi_targets=o_config.qoi_targets)

            o_plot.add_dataframe(
                    color=v['color'],
                    label=v['label'],
                    obj=copy.deepcopy(o_data.df),
                    names=o_data.normalized_names)

    o_plot.make_plot(filename="parallel_plot.png", xlabels=o_data.normalized_error_names, 
                   ylabel="% error", title="Si sw", ylim=(-175, 25), legend_loc="lower right")
