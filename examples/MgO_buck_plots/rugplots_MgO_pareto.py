import copy,os
import numpy as np
import pandas as pd
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
class PyposmatRugplot():

    def __init__(self):
        self.data_directory = None
        self.configuration_fn = None
        self.data_fn = None

        self.qoi_reference_data = None

o_rugplot = PyposmatRugplot()
o_rugplot.data_directory = "../../data/MgO_pareto_data"
o_rugplot.configuration_fn = os.path.join(
        o_rugplot.data_directory,
        'pyposmat.config.in'
)

# -----------------------------------------------------------------------------
# DEFINE WHERE TO FIND ALL THE DATA
# -----------------------------------------------------------------------------
data_directory = "../..//data/MgO_pareto_data"
config_fn = os.path.join(data_directory,'pyposmat.config.in')
datafile_fn = os.path.join(data_directory,'culled_005.out')

# This isn't necessary, but is here if you don't have a pyposmat.config.in file and are going to regenerate it from a script
pyposmat_config_script = None
# pyposmat_config_script = os.path.join(data_directory,"MgO__buck.py")

# -----------------------------------------------------------------------------
# DEFINE WHERE TO PUT ALL THE OUTPUT
# -----------------------------------------------------------------------------
output_directory = "./"
output_plot_fn = os.path.join(
        output_directory,
        'rugplot_MgO_buck.png'
)


qoi_names = ['MgO_NaCl.a0','MgO_NaCl.c11','MgO_NaCl.c12','MgO_NaCl.c44','MgO_NaCl.B','MgO_NaCl.G','MgO_NaCl.fr_a','MgO_NaCl.fr_c','MgO_NaCl.sch','MgO_NaCl.001s']
results_LewisCatlow_raw = [4.21079021691525,307.571810176713,171.13560278936,168.168424864137,216.61433858514434,68.21810369367648,9.679582738895988,9.810034150996216,5.796822583346511,0.06783775861630922]
results_BallGrimes_1_raw = [4.20883001371435,360.106974760923,162.314315016903,160.683383081696,228.24520159824297,98.89632987200999,12.428466047278562,11.87773342645869,7.200953868002216,0.08064679409333339]
results_BallGrimes_2_raw = [4.222448,301.315822058251,150.827961278274,142.471471981922,200.99058153826635,75.2439303899885,10.435732615942925,8.526618652243087,5.509124492308274,0.0692527122209811]

results_LewisCatlow = OrderedDict([])
results_BallGrimes_1 = OrderedDict([])
results_BallGrimes_2 = OrderedDict([])

for i,v in enumerate(qoi_names):
    results_LewisCatlow[v] = results_LewisCatlow_raw[i]
    results_BallGrimes_1[v] = results_BallGrimes_1_raw[i]
    results_BallGrimes_2[v] = results_BallGrimes_2_raw[i]

ref_data = OrderedDict()
ref_data['LC'] = results_LewisCatlow
ref_data['BG1'] = results_BallGrimes_1
ref_data['BG2'] = results_BallGrimes_2

ref_data_colors = OrderedDict()
ref_data_colors['LC'] = "red"
ref_data_colors['BG1'] = 'blue'
ref_data_colors['BG2'] = "green"

o_rugplot.qoi_reference_data = ref_data
o_rugplot.qoi_reference_data_colors = ref_data_colors

def message_out(msg):
    print(msg)

if __name__ == "__main__":
    config=PyposmatConfigurationFile()
  
    # check to see if data directory exists
    if not os.path.isdir(data_directory):
        msg = 'Cannot find data directory, {}\n'.format(data_directory)
        msg += '\t absolute_path:{}\n'.format(os.path.join(data_directory))
        message_out(msg)
        exit()
    else:
        msg = "Found data directory, {}".format(data_directory)
        message_out(msg)

    # read configuration file
    try:
        config.read(filename=config_fn)
    except FileNotFoundError as e:
        msg = "Cannot find pyposmat configuration file:{}".format(config_fn)
        message_out(msg)

        if pyposmat_config_script is not None:
            # run configuration script
            pass
        else:
            msg = "cannot find pyposmat configuration script because pyposmat_config_script variable was not set"
            message_out(msg)
            raise

    # read the data file
    datafile=PyposmatDataFile()
    datafile.read(filename=datafile_fn)
    (nrows,ncols) = datafile.df.shape

    msg = "reading data file....\n"
    msg += "\t{}\n".format(datafile_fn)
    msg += "the data file has...\n"
    msg += "\t{} nrows\n".format(nrows)
    msg += "\t{} ncols\n".format(ncols)
    message_out(msg)

    # set the plot filename
    plot_fn = "rugplots_MgO_buck.png"

    # check to see if we have excluded names
    qoi_excluded_names = []
    qoi_names = [q for q in config.qoi_names if q not in qoi_excluded_names]
    qoi_targets = config.qoi_targets
    
    msg = "original qoi_names is length:{}\n".format(len(config.qoi_names))
    msg += "qoi_excluded_names is length:{}\n".format(len(qoi_excluded_names))
    msg += "qoi_names is length:{}\n".format(len(qoi_names))
    message_out(msg)

    # calculate normed errors
    for iqn,qn in enumerate(qoi_names):
        error_names = "{}.err".format(qn)
        normed_error_names = "{}.nerr".format(qn)
        q = qoi_targets[qn]
        datafile.df[normed_error_names] = datafile.df[qn]/q-1

    normed_error_names = ['{}.nerr'.format(q) for q in qoi_names]

    # calculate distance metric for best potential selection
    datafile.df['d_metric'] = np.sqrt(np.square(datafile.df[normed_error_names]).sum(axis=1))

    # making the plot
    fig,ax = plt.subplots()
    
    # plot data
    for iqn,qn in enumerate(qoi_names):
        x = datafile.df['{}.nerr'.format(qn)]
        y = nrows*[len(qoi_names)-iqn]
        ax.scatter(x,y,
                marker='|',
                #s=10.,
                color='grey')

    # plot reference data
    ref_data_marker_size = 300
    ref_data_marker_type = '|'
    for ref_data_name,ref_data_dict in ref_data.items():
        for iqn,qn in enumerate(qoi_names):
            q = qoi_targets[qn]
            x = ref_data_dict[qn]/q-1
            y = len(qoi_names)-iqn
            ax.scatter(x,y,
                    s=ref_data_marker_size,
                    marker=ref_data_marker_type,
                    color=ref_data_colors[ref_data_name]
            )

    datum_line_size = 0.1
    datum_line_style = '-'
    datum_line_color = 'k'
    plt.axvline(0,
            color=datum_line_color,
            linestyle=datum_line_style,
            linewidth=datum_line_size
    )

    x_label = 'Pct Error Difference'
    ax.set_xlabel(x_label)

    try:
        y_ticks_labels = [config.latex_labels[qn]['name'] for qn in qoi_names]
    except KeyError as e:
        y_ticks_labels = qoi_names

    plt.sca(ax)
    plt.yticks(
            list(range(1,len(qoi_names)+1)),
            list(reversed(y_ticks_labels))
        )
    #vals = ax.get_xticks()
    #ax.set_xticklabels(
    #        ['{:.1%}'.format(int(x*100)) for x in vals]
    #    )

    ax.set_xlim([-1,1])
    ax.set_ylim([0,len(qoi_names)+1])
    ax.legend(
        handles = [
            mpatches.Patch(
                color=ref_data_colors[ref_data_name],
                label=ref_data_name)
            for ref_data_name in ref_data
        ] + [
            mpatches.Patch(
                color='grey',
                label='best 50')
        ]
    )
    fig.savefig(plot_fn)

    exit()
