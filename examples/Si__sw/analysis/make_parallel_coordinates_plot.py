import os, shutil
import copy
from collections import OrderedDict
import pypospack.utils
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.visualization import PyposmatParallelCoordinatesPlot

parallel_plot_config= OrderedDict()
parallel_plot_config['p_3.5_q_0.5'] = OrderedDict()
parallel_plot_config['p_4.0_q_0.0'] = OrderedDict()
parallel_plot_config['p_q_free'   ] = OrderedDict()

parallel_plot_config['p_3.5_q_0.5']['plot_fn'] = os.path.join('parallel_plots','parallel_p_35_q_05.png')
parallel_plot_config['p_4.0_q_0.0']['plot_fn'] = os.path.join('parallel_plots','parallel_p_40_q_00.png')
parallel_plot_config['p_q_free'   ]['plot_fn'] = os.path.join('parallel_plots','parallel_unconstrained.png')


parallel_plot_config['p_3.5_q_0.5']['label'] = 'p=3.5 q=0.5'
parallel_plot_config['p_4.0_q_0.0']['label'] = 'p=4.0 q=0.0'
parallel_plot_config['p_q_free'   ]['label'] = 'unconstrained'

parallel_plot_config['p_3.5_q_0.5']['color'] = 'blue'
parallel_plot_config['p_4.0_q_0.0']['color'] = 'red'
parallel_plot_config['p_q_free'   ]['color'] = 'green'

parallel_plot_config['p_3.5_q_0.5']['data_directory'] = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),'data','Si__sw__data',
            'pareto_optimization_p_3.5_q_0.5')
parallel_plot_config['p_4.0_q_0.0']['data_directory'] = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),'data','Si__sw__data',
            'pareto_optimization_p_4.0_q_0.0')
parallel_plot_config['p_q_free']['data_directory'] = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),'data','Si__sw__data',
            'pareto_optimization_unconstrained')

parallel_plot_config['p_3.5_q_0.5']['data_fn'] = os.path.join(
            parallel_plot_config['p_3.5_q_0.5']['data_directory'],'pyposmat.kde.19.out')
parallel_plot_config['p_4.0_q_0.0']['data_fn'] = os.path.join(
            parallel_plot_config['p_4.0_q_0.0']['data_directory'],'pyposmat.kde.19.out')
parallel_plot_config['p_q_free']['data_fn'] = os.path.join(
            parallel_plot_config['p_q_free']['data_directory'],'pyposmat.kde.19.out')

parallel_plot_config['p_4.0_q_0.0']['config_fn'] = os.path.join(
            parallel_plot_config['p_4.0_q_0.0']['data_directory'],
            'pyposmat.config.in')
parallel_plot_config['p_3.5_q_0.5']['config_fn'] = os.path.join(
            parallel_plot_config['p_3.5_q_0.5']['data_directory'],
            'pyposmat.config.in')
parallel_plot_config['p_q_free']['config_fn'] = os.path.join(
            parallel_plot_config['p_q_free']['data_directory'],
            'pyposmat.config.in')

def determine_ylimits():
    return (-0.5,0.5)

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
    
    figure_fn = os.path.join(plot_dir,'parallel_plot.eps')
    
    o_plot = PyposmatParallelCoordinatesPlot()
    o_plot.create_subplots()
    for k,v in parallel_plot_config.items():
        print(v)
        o_plot.plot(
                config=v['config_fn'],
                data=v['data_fn'],
                label=v['label'],
                color=v['color'],
                nsmallest=20,
                alpha=0.7)

    ref_config_fn = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),'data','Si__sw__data',
            'reference_potentials',
            'pyposmat.config.in')
    ref_data_fn = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),'data','Si__sw__data',
            'reference_potentials',
            'pyposmat.kde.1.out')
    o_plot.plot_reference_potentials(
            config=ref_config_fn,
            data=ref_data_fn,
            linewidth=5)
            
    o_plot.show_figure()
    o_plot.save_figure(filename=plot_fn)
    exit()
    # initialization
    o_plot = PyposmatParallelCoordinatesPlot()

    # add data to plot
    for k,v in parallel_plot_config.items():
        #print(k,v)
        print(80*'-')
        print(k)
        if k == 'args':
            pass
        else:
            o_config = PyposmatConfigurationFile()
            o_config.read(filename=v['config_fn'])
            print(o_config.qoi_targets)
            o_data = PyposmatDataFile()
            o_data.read(filename=v['data_fn'])
            o_data.create_normalized_errors(
                    normalize_type='by_qoi_target',
                    qoi_targets=o_config.qoi_targets)

            o_data.df['score'] = o_data.df[o_config.normalized_error_names].abs().sum(axis=1)
            
            header_fmt = "{:20} {:10} {:10} {:10} {:10} {:10}"
            row_fmt    = "{:20} {:+10.6f}  {:+10.6f}  {:+10.6f}  {:+10.6f}  {:+10.6f}"

            print('errors')
            print(header_fmt.format('qoi_name','target','mean','std','min','max'))
            for qn in o_config.qoi_names:
                en = "{}.err".format(qn)
                q_mean = o_data.df[en].mean()
                q_std = o_data.df[en].std()
                q_min = o_data.df[en].min()
                q_max = o_data.df[en].max()
                q_target = o_config.qoi_targets[qn]
                print(row_fmt.format(en,q_target,q_mean,q_std,q_min,q_max))

            print('normalized errors')
            print(header_fmt.format('qoi_name','target','mean','std','min','max'))
            
            for nen in o_config.normalized_error_names:
                q_mean = o_data.df[nen].mean()
                q_std = o_data.df[nen].std()
                q_min = o_data.df[nen].min()
                q_max = o_data.df[nen].max()
                q_target = 0.
                print(row_fmt.format(nen,q_target,q_mean,q_std,q_min,q_max))

            o_plot.add_dataframe(
                    data=copy.deepcopy(
                        o_data.df.nsmallest(20,'score')[o_config.normalized_error_names]),
                    label=v['label'],
                    color=v['color'],
                    alpha=0.7,
                    names=o_config.normalized_error_names
                    )

            plot_fn = v['plot_fn']
            y_limits = determine_ylimits()
            o_plot.make_plot(
                    filename=v['plot_fn'], 
                    xlabels=o_config.normalized_error_names, 
                    ylabel="% error", 
                    title="Si sw", 
                    ylim= y_limits, 
                    legend_loc="lower right")
