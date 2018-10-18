import copy,os
import numpy as np
import pandas as pd
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.data import PyposmatConfigurationFile

# -----------------------------------------------------------------------------
# Define where to find all the data
# -----------------------------------------------------------------------------

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

class PyposmatParallelCoordinates(object):
    pass

if __name__ == "__main__":
    from pathlib import Path
   
    # define base qoi
    base_qoi = "sim_id"
    # define the data directory, and relevant input files
    data_directory = "../../../data/MgO_pareto_data"
    config_fn = os.path.join(data_directory,'pyposmat.config.in')
    datafile_fn = os.path.join(data_directory,'culled_005.out')

    if Path(data_directory).is_dir():
        print("[OK] data_directory:{}:found".format(data_directory))
    else:
        print("[FAIL] data_directory:{}:not found".format(data_directory))
        exit()

    # define the the output directory, and relevant output files
    output_directory = "./"
    output_plot_fn = os.path.join(output_directory,'rugplot_MgO_buck.eps')


    if Path(data_directory).is_dir():
        print("[OK] data_directory:{}:found".format(data_directory))
    else:
        print("[FAIL] data_directory:{}:not found".format(data_directory))
        exit()
    
    if Path(output_directory).is_dir():
        print("[OK] output_directory:{}:found".format(output_directory))
    else:
        print("[FAIL] output_directory:{}:found".format(output_directory))
        exit()

    # read configuration file
    
    if Path(config_fn).is_file():
        print("[OK] configuration file:{}:found".format(config_fn))
    else:
        print("[FAIL] configuration file:{}:not found".format(config_fn))
        exit()
    
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
