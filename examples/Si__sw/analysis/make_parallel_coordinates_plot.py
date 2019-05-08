import os, shutil
import copy
from collections import OrderedDict
import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.visualization.parallel_plot_new import PyposmatParallelCoordinatesPlot

parallel_plot_config= OrderedDict()
parallel_plot_config['args'] = OrderedDict()
parallel_plot_config['p_3.5_q_0.5'] = OrderedDict()
parallel_plot_config['p_4.0_q_0.0'] = OrderedDict()
parallel_plot_config['p_q_free'] = OrderedDict()

parallel_plot_config['p_3.5_q_0.5']['plot_fn'] = os.path.join('parallel_plots','parallel_p_35_q_05.png')
parallel_plot_config['p_4.0_q_0.0']['plot_fn'] = os.path.join('parallel_plots','parallel_p_40_q_00.png')
parallel_plot_config['p_q_free'   ]['plot_fn'] = os.path.join('parallel_plots','parallel_unconstrained.png')


parallel_plot_config['p_3.5_q_0.5']['label'] = 'p=3.5 q=0.5'
parallel_plot_config['p_4.0_q_0.0']['label'] = 'p=4.0 q=0.0'
parallel_plot_config['p_q_free'   ]['label'] = 'unconstrained'
parallel_plot_config['p_3.5_q_0.5']['color'] = 'blue'
parallel_plot_config['p_3.5_q_0.5']['data_directory'] = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),'data','Si__sw__data',
            'pareto_optimization_p_3.5_q_0.5')
parallel_plot_config['p_3.5_q_0.5']['config_fn'] = os.path.join(
            parallel_plot_config['p_3.5_q_0.5']['data_directory'],
            'pyposmat.config.in')
parallel_plot_config['p_3.5_q_0.5']['data_fn'] = os.path.join(
            parallel_plot_config['p_3.5_q_0.5']['data_directory'],'pyposmat.kde.19.out')

parallel_plot_config['p_4.0_q_0.0']['color'] = 'blue'
parallel_plot_config['p_4.0_q_0.0']['data_directory'] = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),'data','Si__sw__data',
            'pareto_optimization_p_4.0_q_0.0')
parallel_plot_config['p_4.0_q_0.0']['config_fn'] = os.path.join(
            parallel_plot_config['p_4.0_q_0.0']['data_directory'],
            'pyposmat.config.in')
parallel_plot_config['p_4.0_q_0.0']['data_fn'] = os.path.join(
            parallel_plot_config['p_4.0_q_0.0']['data_directory'],'pyposmat.kde.19.out')


parallel_plot_config['p_q_free']['color'] = 'blue'
parallel_plot_config['p_q_free']['data_directory'] = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),'data','Si__sw__data',
            'pareto_optimization_unconstrained')
parallel_plot_config['p_q_free']['config_fn'] = os.path.join(
            parallel_plot_config['p_q_free']['data_directory'],
            'pyposmat.config.in')
parallel_plot_config['p_q_free']['data_fn'] = os.path.join(
            parallel_plot_config['p_q_free']['data_directory'],'pyposmat.kde.19.out')

def determine_ylimits():
    return (-1,1)

def make_latex_table(config,
                     data,
                     qoi_type=None,
                     param_type=None):
    qoi_types = ['by_qoi_target']
    param_type = []

    assert isinstance(config,str) \
           or isinstance(config,PyposmatConfigurationFile)
    assert isinstance(data,str) \
            or isinstance(data,PyposmatDataFile)

    if isinstance(config,str):
        o_config = PyposmatConfigurationFile()
        o_config.read(filename=config)
    elif isinstance(config,PyposmatConfigurationFile):
       o_config = config
    else:
        raise TypeError()

    if isinstance(data,str):
        o_data = PyposmatDataFile()
        o_data.read(filename=data)
    elif isinstance(data,PyposmatDataFile):
        o_data = data
    else:
        raise TypeError()

    if qoi_type == 'by_qoi_target':
        o_data.create_normalized_errors(
                normalize_type='by_qoi_target',
                qoi_targets=o_config.qoi_targets)
        df = o_data.df[o_data.normalized_error_names]


if __name__ == "__main__":
    # reset the plot directory
    plot_dir = 'parallel_plots'
    if os.path.isdir(plot_dir):
        shutil.rmtree(plot_dir)
    os.mkdir(plot_dir)

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
            o_data.create_normalized_errors(
                    normalize_type='by_qoi_target',
                    qoi_targets=o_config.qoi_targets)

            for nen in o_data.normalized_error_names:
                q_mean = o_data.df[nen].mean()
                q_std = o_data.df[nen].std()
                q_min = o_data.df[nen].min()
                q_max = o_data.df[nen].max()

                print(nen,q_mean,q_std,q_min,q_max)
            #if ymin is None:
            #    ymin = o_data.df[o_data.normalized_error_names].min()
            #else:
            #    ymin = min(ymin,o_data.df[o_data.normalized_error_names].min())

            #if ymax is None:
            #    ymax = o_data.df[o_data.normalized_error_names].max()
            #else:
            #    ymax = max(ymax,o_data.df[o_data.normalized_error_names].max())

            o_plot.add_dataframe(
                color=v['color'],
                label=v['label'],
                obj=copy.deepcopy(o_data.df),
                names=o_data.normalized_names)

            plot_fn = v['plot_fn']
            y_limits = determine_ylimits()
            o_plot.make_plot(
                    filename=plot_fn, 
                    xlabels=o_data.normalized_error_names, 
                    ylabel="% error", 
                    title="Si sw", 
                    ylim= y_limits, 
                    legend_loc="lower right")
