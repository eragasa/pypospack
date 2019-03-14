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
        data_fn = os.path.join(data_directory,'pyposmat.kde.{}.out'.format(i+1))
        plot_fn = os.path.join(data_directory,"rugplot_{}.png".format(i))

        make_rug_plot(config_fn=config_fn,
                      data_fn=data_fn,
                      plot_fn=plot_fn)
    if False:
        from pypospack.pareto import pareto

        df = copy.deepcopy(datafile.df)
        nr,nc = df.shape
        _nsimulations = OrderedDict()
        _nsimulations['start'] = nr
        abs_error_names = ["{}.abserr".format(q) for q in datafile.qoi_names]
        for q in datafile.qoi_names:
            qe = "{}.err".format(q)
            qne = "{}.abserr".format(q)
            df[qne] = df[qe].abs()
        names = list(df.columns.values)
        abs_err_idx = [names.index(n) for n in abs_error_names]
        pareto_idx = pareto(df[abs_error_names].values.tolist())
        datafile.df = df.loc[pareto_idx,datafile.names]
        datafile.write("results.pareto.out")

        #pareto_idx = pareto_bruteforce(df[abs_error_names].values.tolist())

            #print(pareto_set)
        if make_rugplots:
            rugplots_fn = "rugplot.png"
            datafile = PyposmatDataFile()
            datafile.read("results.pareto.out")
            datafile.qoi_references = OrderedDict()
            datafile.qoi_references['TARGET'] = copy.deepcopy(qoi_targets)
            datafile.score_by_d_metric(scaling_factors='TARGET')
            datafile.subselect_by_score(
                    score_name='d_metric',
                    n=_n_potentials)
            subselect_fn = datafile.write_subselect()

            datafile=PyposmatDataFile()
            datafile.read(filename=_data_fn)
            print(datafile.df[datafile.error_names])
            make_rug_plot(o_config=config,
                          o_data=datafile,
                          fn=rugplots_fn)
