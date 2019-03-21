import os,copy,argparse
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataAnalyzer

def get_qoi_targets(o_config):
    print(type(o_config))
    assert type(o_config) is PyposmatConfigurationFile
    return OrderedDict([
        (k,v['target']) for k,v in o_config.qois.items()]
    )

def show_qoi_targets(config_fn,
                     data_fn):

    o_config = PyposmatConfigurationFile()
    o_config.read(filename=config_fn)

    o_data = PyposmatDataFile()
    o_data.read(filename=data_fn)

    for qoi_name, qoi_target in o_config.qoi_targets.items():
        try:
            qoi_avg = o_data.df[qoi_name].mean()
        except KeyError as e:
            qoi_avg = 'no value'
        s = "{:20} {:10} {:10}".format(qoi_name,qoi_target,qoi_avg)
        print(s)


import matplotlib.pyplot as plt
def make_rug_plot(config_fn,
                  data_fn,
                  ax=None,
                  plot_fn='rugplot.png'):

    o_config = PyposmatConfigurationFile()
    o_config.read(filename=config_fn)

    o_data = PyposmatDataFile()
    o_data.read(filename=data_fn)

    qoi_targets = o_config.qoi_targets
    #qoi_targets = get_qoi_targets(o_config)
    error_names = o_data.error_names
    qoi_names = o_data.qoi_names

    # create normalized error
    df = copy.deepcopy(o_data.df[error_names])
    for qn in qoi_names:
        en = "{}.err".format(qn)
        nen = "{}.nerr".format(qn)
        q = qoi_targets[qn]
        df[nen]=o_data.df[en]/q-q

    (_nrows,_ncols) = o_data.df.shape

    if ax is None:
        fig, ax = plt.subplots(nrows=1,ncols=1)

    for iq,qn in enumerate(qoi_names):
        _yloc = [iq+1]
        ax.scatter(
            df["{}.nerr".format(qn)],
            _nrows*[iq+1],
            marker='|',
            s=100.,
            color='k'
        )

    plt.sca(ax)
    plt.yticks(range(len(qoi_names)+1),['']+qoi_names)
    fig.savefig(plot_fn)



if __name__ == "__main__":
    data_directory = 'data'
    plot_directory = 'plots'
    n_iterations = 10
    n_potentials = 30
   
    if not os.path.isdir(plot_directory):
        os.mkdir(plot_directory)
    for i in range(n_iterations):
        config_fn = os.path.join(data_directory,'pyposmat.config.in')
        results_data_fn = os.path.join(data_directory,'pyposmat.results.{}.out'.format(i))
        kde_data_fn = os.path.join(data_directory,'pyposmat.kde.{}.out'.format(i+1))
        plot_fn = os.path.join(plot_directory,"rugplot_{}.png".format(i))

        print(80*'=')
        print("{:^80}".format('ITERATION {}'.format(i)))

        results_data = PyposmatDataFile()
        results_data.read(filename=results_data_fn)
        results_n_rows,results_n_cols = results_data.df.shape

        kde_data = PyposmatDataFile()
        kde_data.read(filename=kde_data_fn)
        kde_n_rows,kde_n_cols=kde_data.df.shape
        print('total_number_of_candidates:{}'.format(results_n_rows))
        print('remaining_number_of_candiddates:{}'.format(kde_n_rows))

        show_qoi_targets(config_fn=config_fn,data_fn=kde_data_fn)

        make_rug_plot(config_fn=config_fn,
                      data_fn=kde_data_fn,
                      plot_fn=plot_fn)
