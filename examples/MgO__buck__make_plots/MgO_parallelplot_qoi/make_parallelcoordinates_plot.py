import copy,os
import numpy as np
import pandas as pd
from collections import OrderedDict
from pypospack.pyposmat.visualization.parallel_plot_qoi import PyposmatQoiParallelCoordinatesPlot

pypospack_root_dir = [v for v in os.environ['PYTHONPATH'].split(':') if v.endswith('pypospack')][0]
# -----------------------------------------------------------------------------
# DEFINE WHERE TO FIND ALL THE DATA
# -----------------------------------------------------------------------------
datafile_fn = os.path.join(pypospack_root_dir,'data/MgO_pareto_data/qoiplus_005.out')
config_fn = os.path.join(pypospack_root_dir,'examples/MgO__buck__add_additional_qoi/data/pyposmat.config.in')

# -----------------------------------------------------------------------------
# DEFINE WHERE TO PUT ALL THE OUTPUT
# -----------------------------------------------------------------------------
output_directory = "./"
output_plot_fn = os.path.join(
        output_directory,
        'qoi_parallelplot_MgO_buck.png'
)

if __name__ == "__main__":
    o_plot = PyposmatQoiParallelCoordinatesPlot()
    o_plot.read_configuration(filename=config_fn)
    o_plot.read_datafile(filename=datafile_fn)
    o_plot.plot_legend_location = 'best'
    o_plot.make_plot(
            filename=output_plot_fn,
            include_qois=True,
            include_qois_v=True,
            qoi_excluded_names=None)

    exit()
    output_plot_fn = os.path.join(output_directory,'rugplot_MgO_buck.eps')


    
    config=PyposmatConfigurationFile()
    config.read(filename=config_fn)

    # read the associated datafile
    
    if Path(datafile_fn).is_file():
        print("[OK] data file:{}:found".format(datafile_fn))
    else:
        print("[FAIL] data file:{}:not found".format(datafile_fn))
        exit()
    
    datafile=PyposmatDataFile()
    datafile.read(filename=datafile_fn)

    # output d
    output_directory = "./"
    plot_fn = os.path.join(output_directory,'parallelcoordinates_fs.png')

    excluded_qoi_names = []
    qoi_names = [q for q in config.qoi_names if q not in excluded_qoi_names]

    if excluded_qoi_names is []:
        print('no excluded quantities of interest')

    print(80*'=')
    print('QOI_NAMES')
    print(80*'=')
    for qn in qoi_names:
        print('\t{}'.format(qn))
    print('qoi_names is length:{}'.format(len(qoi_names)))
    error_names = ["{}.err".format(q) for q in qoi_names]
    normed_error_names = ["{}.nerr".format(q) for q in qoi_names]
    qoi_targets = config.qoi_targets

    # calculate normalized error for sampled data
    for iqn,qn in enumerate(qoi_names):
        en = "{}.err".format(qn)
        nen = "{}.nerr".format(qn)
        q = qoi_targets[qn]
        datafile.df[nen] = datafile.df[qn]/q-1
    (nrows,ncols) = datafile.df.shape

    normederr_names = ['{}.nerr'.format(q) for q in qoi_names]
    datafile.df['d_metric'] = np.sqrt(np.square(datafile.df[normederr_names]).sum(axis=1))
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    from pandas.plotting import parallel_coordinates
    fig, ax = plt.subplots()


    reference_df = pd.DataFrame(list(ref_data.values()))
    for iqn,qn in enumerate(qoi_names):
        en = "{}.err".format(qn)
        nen = "{}.nerr".format(qn)
        q = qoi_targets[qn]
        reference_df[nen] = reference_df[qn]/q-1
    reference_df['sim_id'] = list(ref_data.keys())
    #reference_df.set_index('sim_id') 

    data_df = copy.deepcopy(datafile.df)
    #data_df.set_index('sim_id')
    print(data_df[['sim_id'] + normed_error_names].columns)
    
    subselect_df = datafile.df.nsmallest(30,'d_metric')
    #subselect_df.set_index('sim_id')
    is_plot_all_data = False
    if is_plot_all_data:
        print("data_df:{}".format(data_df[normed_error_names].shape))

        column_names = ['sim_id'] + normed_error_names
        parallel_coordinates(
                column_names,
                base_qoi,
                color='grey'
                )

    subselect_color = 'grey'
    is_plot_subselect = True
    if is_plot_subselect:
        print('subselect_df:{}'.format(
            subselect_df[['sim_id'] + normed_error_names].shape)
            )
        column_names = ['sim_id'] + normed_error_names
        parallel_coordinates(
                subselect_df[column_names],
                base_qoi,
                color=subselect_color)

    # plot reference data

    column_names = ['sim_id'] + normed_error_names
    parallel_coordinates(
            reference_df[column_names],
            base_qoi,
            )


    plt.gca().legend_.remove()

    plt.xticks(rotation=90)
    #fig.savefig(plot_fn)
    fig.tight_layout()
    plt.show()
    exit()

    (nr,nc)=df.shape
    print("We have {} potentials...".format(nr))
    for iqn,qn in enumerate(qoi_names):
        nen = '{}.nerr'.format(qn)
        x = df[nen]
        y = nr*[len(qoi_names)-iqn]
        ax.scatter(x,y,
            marker='|',
            #s=10.,
            color='grey')

    for ref_data_name,ref_data_dict in ref_data.items():
        for iqn,qn in enumerate(qoi_names):
            q = qoi_targets[qn]
            x = ref_data_dict[qn]/q-1
            y = len(qoi_names)-iqn
            ax.scatter(
                    x,
                    y,
                    s=300,
                    marker='|',
                    color=ref_data_colors[ref_data_name]
            )

    plt.axvline(0,color='k',linestyle='-',linewidth=.1)
    ax.set_xlabel('Pct Error Difference')
    yticks_labels = [config.latex_labels[qn]['name'] for qn in qoi_names]
    print("length of yticks_labels:{}".format(len(yticks_labels)))
    plt.sca(ax)
    plt.yticks(
            list(range(1,len(qoi_names)+1)),
            list(reversed(yticks_labels))
        )
    ax.set_xticks(rotation=90)
    plt.show()
    #vals = ax.get_xticks()
    #ax.set_xticklabels(
    #        ['{:.1%}'.format(int(x*100)) for x in vals]
