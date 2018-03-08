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
def make_rug_plot(o_config,o_data,ax=None,fn='rugplot.png'):
    qoi_targets = get_qoi_targets(o_config)
    error_names = o_data.error_names
    qoi_names = o_data.qoi_names

    # create normalized errors
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
            _nrows*[i_error+1],
            marker='|',
            s=100.,
            color='k'
        )
    plt.sca(ax)
    plt.yticks(range(len(qoi_names)+1),['']+qoi_names)
    fig.savefig(fn)

if __name__ == "__main__":
    #parser = argparse.ArgumentParser()
    #parser.add_argument("--src_dir",
    #        action='store',
    #        dest='data_dir',
    #        type=str,
    #        help="location of the data directory")
    #parse_args = parser.parse_args()

    _n_potentials = 30
    _data_fn = "data__Ni__eam__born_exp_bjs_00\pyposmat.kde.10.out"
    _config_fn = "data__Ni__eam__born_exp_bjs_00\pyposmat.config.in"

    print(80*'-')
    print("reading the configuration file {}...".format(_config_fn))
    config=PyposmatConfigurationFile()
    config.read(filename=_config_fn)
    qoi_targets=get_qoi_targets(config)
    print("reading the data file {}...".format(_data_fn))
    datafile=PyposmatDataFile()
    datafile.read(filename=_config_fn)
    datafile.qoi_references = OrderedDict()
    datafile.qoi_references['TARGET'] = copy.deepcopy(qoi_targets)
    datafile.score_by_d_metric(scaling_factors='TARGET')
    datafile.subselect_by_score(
            score_name='d_metric',
            n=_n_potentials)
    subselect_fn = datafile.write_subselect()

    datafile=PyposmatDataFile()
    datafile.read(filename=subselect_fn)

    make_rug_plot(o_config=config,o_data=datafile)
